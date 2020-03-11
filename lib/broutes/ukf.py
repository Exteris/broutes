"""
Perform Unscented Kalman Filtering on GPS data to get a smoother boat velocity
and possibly retrieve the stroke rate.

The first model is simple, with state variables
lat: latitude (we do not correct for wraparound, good luck near the poles)
v_north: velocity in the northern direction
a_north: acceleration in the northern direction
lon: longitude (we do not correct for wraparound, good luck near the meridian)
v_east: velocity in the eastern direction
a_east: acceleration in the eastern direction

Additionally we would like to model rowing motion back and forth, by adding
a sinusoidal velocity perturbation. This is described by three variables, the
l: stroke reach (0-2 meters)
phi: rowing phase (radians, always increasing, > 0)
R: rowing frequency (strokes per second, 0..2/3)
"""
import sys
import json
from datetime import datetime
import numpy as np
from numpy.linalg import norm

import matplotlib.pyplot as plt
plt.rc('text', usetex=False)

from filterpy.common import Q_discrete_white_noise
from filterpy.kalman import MerweScaledSigmaPoints
from filterpy.kalman import UnscentedKalmanFilter as UKF

def read_data(filename):
    """Returns an array of time offsets (1d) and gps coordinates (2d)"""
    with open(filename, 'r') as file:
        gps_data = json.load(file)
    start_time = parse_time(gps_data['points'][0]['time'])
    times = np.asarray([(parse_time(p['time']) - start_time).total_seconds() for p in gps_data['points'] if has_gps(p)])
    locations = np.asarray([(p['lat'], p['lon']) for p in gps_data['points'] if has_gps(p)])

    # Convert the locations to 'meters' by multiplying with the right constant
    return times, locations * 111111

def parse_time(s):
    """Parse the time as printed by Ruby"""
    return datetime.strptime(s, '%Y-%m-%d %H:%M:%S +0000')

def has_gps(point):
    """Is there a latitude and longitude in a point (dict?)"""
    return 'lat' in point and 'lon' in point

def f_motion(x, dt):
    """State transition function for simple directed motion"""

    # Constrain accelerations and velocities
    if False:
        vnorm = norm(x[[1,4]], ord=2)
        anorm = norm(x[[2,5]], ord=2)
        if vnorm > 7: # m/s
            x[1] = x[1] / vnorm * 7
            x[4] = x[4] / vnorm * 7

        if anorm > 2: # m/s^2 (free fall = 9.81), so this is ~0.2 g
            x[2] = x[2] / anorm * 2
            x[5] = x[5] / anorm * 2

    F = np.array([[1, dt, 0.5*dt**2, 0, 0, 0],
                  [0, 1, dt,         0, 0, 0],
                  [0, 0,  1,         0, 0, 0],
                  [0, 0,  1, 1, dt, 0.5*dt**2],
                  [0, 0,  0, 0,  1, dt],
                  [0, 0,  0, 0,  0,  1]])
    return np.dot(F, x)

def f_rowing_motion(x_in, dt):
    """State transition function for rowing motion for wrist-worn gps receivers.

    We have additional variables phi and R here, with a model for an additional velocity
    perturbation (outside of the acceleration term) in the current direction of motion.

    X is then interpreted as the position of the gps receiver on the athlete wrist,
    and v is interpreted as the average boat velocity, and a the boat acceleration,
    minus any component caused by the stroke motion.

    Extra equations:
    x_o = x_in + ... + v/|v| * l * cos(phi)
    phi = 2 * pi * R * dt
    l = l
    R = R

    constraints imposed:
    0 < |a| < 0.5 m/s^2 (~ 0.05 g)
    0 < |v| < 7 m/s (~ 25 km/h)
    0 < R < 2/3 Hz (40 strokes per minute)
    0 < phi
    0 < l < 2 m


    """
    x, v_x, a_x, y, v_y, a_y, phi, R = x_in

    l = 1.5

    v_norm = max(np.sqrt(v_x**2 + v_y**2), 1e-8)
    if v_norm > 7 and False: # m/s
        v_x = v_x / v_norm * 7
        v_y = v_y / v_norm * 7

    a_norm = max(np.sqrt(a_x**2 + a_y**2), 1e-8)
    if a_norm > 0.5 and False:
        a_x = a_x / a_norm * 0.5
        a_y = a_y / a_norm * 0.5


    #l = np.clip(l, 0, 2)
    #R = np.clip(R, 1e-8, 2./3.)

    phi_o = phi + 2 * np.pi * R * dt

    x_o   = x + v_x * dt + 0.5 * a_x * dt**2 \
          + v_x/v_norm * l * (np.cos(phi_o) - np.cos(phi))
    v_x_o = v_x + a_x * dt
    a_x_o = a_x
    y_o   = y + v_y * dt + 0.5 * a_y * dt**2 \
          + v_y/v_norm * l * (np.cos(phi_o) - np.cos(phi))
    v_y_o = v_y + a_y * dt
    a_y_o = a_y
    R_o = R

    return np.asarray([x_o, v_x_o, a_x_o, y_o, v_y_o, a_y_o, phi_o, R_o])

def h_gps(x):
    """Just return the GPS coordinates as measurement function"""
    return x[[0, 3]]

def batch_ukf_motion(times, zs):
    """Calculate a batch process UKF for the simple motion case"""
    sel = np.s_[:]

    points = MerweScaledSigmaPoints(n=6, alpha=.1, beta=2., kappa=-3)
    dt = times[1] - times[0] # assumed constant!
    ukf = UKF(dim_x=6, dim_z=2, fx=f_motion, hx=h_gps, dt=dt, points=points)

    ukf.x = np.array([zs[0, 0], 0., 0., zs[0, 1], 0., 0.])
    # assume a 3 meter gps accuracy
    # no cross-correlation in GPS errors
    gps_std = 3 # meters
    ukf.R = np.diag([gps_std**2, gps_std**2])
    ukf.Q[0:3, 0:3] = Q_discrete_white_noise(3, dt=dt, var=3**2) # meters
    ukf.Q[3:6, 3:6] = Q_discrete_white_noise(3, dt=dt, var=3**2) # meters
    uxs = []
    Xs, Ps = ukf.batch_filter(zs[sel])
    # Ms, P, K = ukf.rts_smoother(Xs, Ps)
    uxs = np.array(Xs)

    sel = np.s_[:]
    plt.plot(uxs[sel, 0], uxs[sel, 3])
    plt.plot(zs[sel, 0], zs[sel, 1], 'o', markersize=2)
    plt.show()

    plt.plot(times[sel], norm(uxs[sel, [1, 4]], ord=2, axis=1))
    distances = l2_distance(zs[:-1, 0], zs[1:, 0], zs[:-1, 1], zs[1:, 1])
    plt.plot(times[1:][sel], distances[sel]/dt, 'o', markersize=2)
    plt.show()

def batch_ukf_rowing_motion(times, zs):
    """Calculate a batch process UKF for the rowing motion case"""
    sel = np.s_[3100:3400]

    points = MerweScaledSigmaPoints(n=8, alpha=.1, beta=2., kappa=-5)
    dt = times[1] - times[0] # assumed constant!
    ukf = UKF(dim_x=8, dim_z=2, fx=f_rowing_motion, hx=h_gps, dt=dt, points=points)

    ukf.x = np.array([zs[sel, 0][0], 0., 0., zs[sel, 1][0], 0., 0., 0., 0.])
    # assume a 3 meter gps accuracy
    # no cross-correlation in GPS errors
    gps_std = 4 # meters
    ukf.R = np.diag([gps_std**2, gps_std**2])
    ukf.Q[0:3, 0:3] = Q_discrete_white_noise(3, dt=dt, var=1**2) # meters
    ukf.Q[3:6, 3:6] = Q_discrete_white_noise(3, dt=dt, var=1**2) # meters
    ukf.Q[6:8, 6:8] = Q_discrete_white_noise(2, dt=dt, var=0.5**2) # radians
    # Add a covariance from phase to position
    ukf.Q[0,6] = 0.5**2
    ukf.Q[6,0] = 0.5**2

    Xs, Ps = ukf.batch_filter(zs[sel])
    Ms, P, K = ukf.rts_smoother(Xs, Ps)
    uxs = np.array(Ms)

    fig, ax = plt.subplots(8, 1, sharex=True, gridspec_kw={'hspace': 0})
    ax[-1].set_xlabel('Time [s]')

    # x error
    ax[0].plot(times[sel], uxs[:, 0] - zs[sel, 0])
    ax[0].set_ylabel('x - x_gps [m]')
    ax[0].axhline(0)

    # v_x
    ax[1].plot(times[sel], uxs[:, 1])
    ax[1].plot(times[sel][1:], np.diff(zs[sel, 0]), 'o', markersize=2)
    ax[1].set_ylabel('v_x [m/s]')

    # a_x
    ax[2].plot(times[sel], uxs[:, 2])
    ax[2].plot(times[sel][2:], np.diff(zs[sel, 0], n=2), 'o', markersize=2)
    ax[2].set_ylabel('a_x [m/s^2]')

    # y error
    ax[3].plot(times[sel], uxs[:, 3] - zs[sel, 1])
    ax[3].set_ylabel('y - y_gps [m]')
    ax[3].axhline(0)

    # v_y
    ax[4].plot(times[sel], uxs[:, 4])
    ax[4].plot(times[sel][1:], np.diff(zs[sel, 1]), 'o', markersize=2)
    ax[4].set_ylabel('v_y [m/s]')

    # a_y
    ax[5].plot(times[sel], uxs[:, 5])
    ax[5].plot(times[sel][2:], np.diff(zs[sel, 1], n=2), 'o', markersize=2)
    ax[5].set_ylabel('a_y [m/s^2]')

    # l
    #ax[4].plot(times[sel], uxs[:, 6])
    #ax[4].set_ylabel('l [m]')

    # phi
    ax[6].plot(times[sel], uxs[:, 6])
    ax[6].set_ylabel('phi [rad]')
    # R
    ax[7].plot(times[sel], uxs[:, 7] * 60)
    ax[7].set_ylabel('Rate [spm]')

    plt.show()

def l2_distance(x1, x2, y1, y2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)

def haversine_distance(lat1, lat2, lon1, lon2):
    """Calculate the haversine distance between two coordinates"""
    RAD_PER_DEG = 0.017453293 #  PI/180
    EARTH_RADIUS = 6_371_000 # m

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    dlon_rad = dlon * RAD_PER_DEG
    dlat_rad = dlat * RAD_PER_DEG

    lat1_rad = lat1 * RAD_PER_DEG
    lon1_rad = lon1 * RAD_PER_DEG

    lat2_rad = lat2 * RAD_PER_DEG
    lon2_rad = lon2 * RAD_PER_DEG

    a = (np.sin(dlat_rad/2))**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * (np.sin(dlon_rad/2))**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))

    return EARTH_RADIUS * c

if __name__ == '__main__':
    times, points = read_data(sys.argv[1])
    batch_ukf_rowing_motion(times, points)
