module Broutes
  class GeoRoute

    attr_reader :start_point, :end_point, :started_at, :ended_at, :total_time
    attr_writer :total_distance, :started_at
    attr_accessor :total_time

    def points
      get_points.to_enum
    end

    def laps
      get_laps.to_enum
    end

    def initialize(args={})
      args.each_pair do |key, value|
        if key.to_sym == :points
          value.each { |p| add_point(p) }
        elsif key.to_sym == :laps
          value.each { |l| add_lap(l) }
        else
          send("#{key}=", value) if respond_to?("#{key}=")
        end
      end
    end

    def self.from_hash(h)
      route = GeoRoute.new h
    end

    def to_hash
      h = {
        'total_distance' => total_distance,
        'total_time' => total_time,
        'total_ascent' => total_ascent,
        'total_descent' => total_ascent,
        'points' => points.collect { |p| p.to_hash },
        'laps' => laps.collect { |l| l.to_hash }
      }
      h['start_point'] = start_point.to_hash if start_point
      h['end_point'] = end_point.to_hash if end_point
      h['started_at'] = @started_at if @started_at
      h
    end

    def add_point(args)
      point = GeoPoint.new(args)
<<<<<<< Updated upstream
      if @start_point
        if distance = Maths.haversine_distance(@end_point, point)
          @total_distance += distance
        end

        @total_time = point.time - @start_point.time if point.time
      else
        @start_point = point
        @end_point = point
        @total_distance = 0
      end

      point.distance = @total_distance
      process_elevation_delta(@end_point, point)
=======
>>>>>>> Stashed changes

      if point.has_location?
        @end_point = point
      end
      get_points << point
    end

    def add_lap(args)
      get_laps << Lap.new(args)
    end

    def process_elevation_delta(last_point, new_point)
      if last_point && new_point && last_point.elevation && new_point.elevation
        delta = new_point.elevation - last_point.elevation
        @_total_ascent = self.total_ascent + (delta > 0 ? delta : 0)
        @_total_descent = self.total_descent - (delta < 0 ? delta : 0)
      end
    end

    def started_at
      return @started_at if @started_at
      @start_point.time if @start_point
    end

    def ended_at
      @end_point.time if @end_point
    end

    def total_ascent
      @_total_ascent ||= 0
    end

    def total_descent
      @_total_descent ||= 0
    end

    # Public : Total distance measured between points in whole metres
    #
    # Returns Float distance in m
    def total_distance
      @total_distance.round if @total_distance
    end

    # Public : Measure of how hilly the route is. Measured as total ascent (m) / distance (km)
    #
    # Returns Float measure
    def hilliness
      (total_distance > 0) ? (total_ascent * 1000 / total_distance) : 0
    end

    # Public: Get average heart_rate for whole GeoRoute.
    #
    # Examples
    #   @route.average_heart_rate
    #   # => 12
    #
    # Returns Integer average, or 0 if no heart_rate on points.
    def average_heart_rate
      points = @_points
      points.map { |p| p.heart_rate || 0 }.inject { |sum, hr| sum + hr } / points.count
    end

    # Public: Get average power for whole GeoRoute.
    #
    # Examples
    #   @route.average_power
    #   # => 250
    #
    # Returns Float average, or 0 if no power on points.
    def average_power
      points = @_points
      points.map { |p| p.power || 0 }.inject { |sum, p| sum + p } / points.count
    end

    # Public: Get maximum speed for whole GeoRoute.
    #
    # Examples
    #   @route.maximum_speed
    #   # => 5.50
    #
    # Returns Float maximum, or 0.0 if no speed on points.
    def maximum_speed
      points = @_points
      points.map { |p| p.speed }.compact.max || 0.0
    end

    # Public: Get minimum speed for whole GeoRoute.
    #
    # Examples
    #   @route.minimum_speed
    #   # => 1.50
    #
    # Returns Float minimum, or 0.0 if no speed on points.
    def minimum_speed
      points = @_points
      points.map { |p| p.speed }.compact.min || 0.0
    end

    # Public: Get average speed for whole GeoRoute.
    #
    # Examples
    #   @route.average_speed
    #   # => 2.50
    #
    # Returns Float average, or 0.0 if no speed on points.
    def average_speed
      points = @_points
      points.map { |p| p.speed || 0.0 }.inject { |sum, p| sum + p } / points.count
    end

    # Public: Get maximum heart rate for whole GeoRoute.
    #
    # Examples
    #   @route.maximum_heart_rate
    #   # => 180
    #
    # Returns Integer maximum, or 0 if no heart rate on points.
    def maximum_heart_rate
      points = @_points
      points.map { |p| p.heart_rate }.compact.max || 0
    end

    # Public: Get minimum heart rate for whole GeoRoute.
    #
    # Examples
    #   @route.minimum_heart_rate
    #   # => 100
    #
    # Returns Integer minimum, or 0 if no heart rate on points.
    def minimum_heart_rate
      points = @_points
      points.map { |p| p.heart_rate }.compact.min || 0
    end

    # Public: Get maximum elevation for whole GeoRoute.
    #
    # Examples
    #   @route.maximum_elevation
    #   # => 1000
    #
    # Returns Integer maximum, or 0 if no elevation on points.
    def maximum_elevation
      points = @_points
      points.map { |p| p.elevation }.compact.max || 0
    end

    # Public: Get minimum elevation for whole GeoRoute.
    #
    # Examples
    #   @route.minimum_elevation
    #   # => 10
    #
    # Returns Integer minimum, or 0 if no elevation on points.
    def minimum_elevation
      points = @_points
      points.map { |p| p.elevation }.compact.min || 0
    end

    # Public: Get total calories for whole GeoRoute.
    #
    # Examples
    #   @route.total_calories
    #   # => 10
    #
    # Returns Integer calories, 0 if no calories on laps or no laps on the route.
    def total_calories
      laps.map { |l| l.calories }.inject { |sum, l| sum + l } || 0
    end

    # Public: Clear outlying points (gps only for now)
    # This uses the median and the median absolute deviation to remove points
    # more than 5 MAD's away from the median point.
    # Removed points are set to nil
    # This is the same algorithm as used in GPXsee (https://github.com/tumic0/GPXSee/blob/master/src/data/track.cpp)
    # This does nothing if acceleration is all empty (i.e. if interpolate_speed_acceleration has not been called)
    # It removes the acceleration, speed and gps coordinates of that point
    #
    # Examples
    #   @route.clear_outliers!
    #
    # Returns nil
    def clear_outliers!
      acceleration = get_points.map(&:acceleration).compact
      m = Maths.median(acceleration)
      mad = Maths.MAD(acceleration, m)
      get_points.map!.with_index do |p, i|
        if p.acceleration.nil? or ((0.6745 * (p.acceleration - m)) / mad).abs > 5
          p.acceleration = nil
          p.speed = nil
          p.lat = nil
          p.lon = nil
          p.elevation = nil
          p.distance = nil
        end

        p
      end
      nil
    end

    # Public: Calculate the speed and acceleration from the GPS data
    # This deletes any existing speed and acceleration data!
    #
    # Examples
    #   @route.calculate_speed_acceleration!
    #
    # Returns nil
    def calculate_speed_acceleration!(overwrite_speed=true)
      i_last_proper = get_points.index(&:has_location?)
      if i_last_proper.nil? # no location data
        return nil
      end
      (0..(i_last_proper-1)).each do |i|
        get_points[i].speed = nil
        get_points[i].acceleration = nil
      end
      get_points[i_last_proper].speed = 0
      get_points[i_last_proper].acceleration = 0

      (i_last_proper+1..get_points.length-1).each do |i|
        p = get_points[i]
        p_prev = get_points[i_last_proper]
        if not p_prev.has_location?
          puts 'ERROR', i
        end
        if ds = Maths.haversine_distance(p, p_prev)
          dt = p.time.to_f - p_prev.time.to_f

          if dt < 1e-3
            p.speed = p_prev.speed if overwrite_speed or p.speed.nil?
            p.acceleration = p_prev.acceleration
          else
            p.speed = ds/dt if overwrite_speed or p.speed.nil?
            dv = p.speed - p_prev.speed
            p.acceleration = dv/dt
          end
        else
          p.speed = p_prev.speed if overwrite_speed or p.speed.nil?
          p.acceleration = nil
        end

        if p.has_location?
          i_last_proper = i
        end

        p
      end
      nil
    end

    # Public: Recalculates the distance
    #
    # Examples
    #   @route.recalculate_distance!
    #
    # Returns nil
    def recalculate_distance!
      total_distance = 0
      @end_point = get_points[0]
      @start_point = get_points[0]
      @_total_distance = 0
      @_total_ascent = 0
      @_total_descent = 0

      get_points.map! do |point|
        process_elevation_delta(@end_point, point)
        if distance = Maths.haversine_distance(@end_point, point)
          @_total_distance += distance
          @end_point = point
        end

        @total_time = point.time - @start_point.time if point.time

        point
      end
      nil
    end

    private

    def get_points
      @_points ||= []
    end

    def get_laps
      @_laps ||= []
    end
  end
end
