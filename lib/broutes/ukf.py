"""
Perform Unscented Kalman Filtering on GPS data to get a smoother boat velocity
and possibly retrieve the stroke rate.

The model is simple, with state variables
lat: latitude (we do not correct for wraparound, good luck near the poles)
v_north: velocity in the northern direction
a_north: acceleration in the northern direction
lon: longitude (we do not correct for wraparound, good luck near the meridian)
v_east: velocity in the eastern direction
a_east: acceleration in the eastern direction
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
            x[2] = x[2] / anorm * 7
            x[5] = x[5] / anorm * 7

    # There are 111111 meters per degree latitude and approximately per degree longitude
    F = np.array([[1, dt, 0.5*dt**2, 0, 0, 0],
                  [0, 1, dt,         0, 0, 0],
                  [0, 0,  1,         0, 0, 0],
                  [0, 0,  1, 1, dt, 0.5*dt**2],
                  [0, 0,  0, 0,  1, dt],
                  [0, 0,  0, 0,  0,  1]])
    return np.dot(F, x)

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
    batch_ukf_motion(times, points)
