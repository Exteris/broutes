# frozen_string_literal: true

module Broutes
  require 'lowess'

  class GeoRoute
    attr_reader :start_point, :end_point, :started_at, :ended_at
    attr_writer :total_distance, :started_at
    attr_accessor :total_time

    def points
      get_points.to_enum
    end

    def laps
      get_laps.to_enum
    end

    def initialize(args = {})
      args.each_pair do |key, value|
        if key.to_sym == :points
          value.each { |p| add_point(p) }
        elsif key.to_sym == :laps
          value.each { |l| add_lap(l) }
        else
          send("#{key}=", value) if respond_to?("#{key}=")
        end
        {
          'total_distance' => total_distance,
          'total_time' => total_time,
          'total_ascent' => total_ascent,
          'total_descent' => total_ascent,
          'points' => points.collect(&:to_hash),
          'laps' => laps.collect(&:to_hash)
        }
      end
    end

    def self.from_hash(h)
      GeoRoute.new h
    end

    def to_hash
      h = {
        'total_distance' => total_distance,
        'total_time' => total_time,
        'total_ascent' => total_ascent,
        'total_descent' => total_ascent,
        'points' => points.collect(&:to_hash),
        'laps' => laps.collect(&:to_hash)
      }
      h['start_point'] = start_point.to_hash if start_point
      h['end_point'] = end_point.to_hash if end_point
      h['started_at'] = @started_at if @started_at
      h
    end

    def summary
      {
        'total_distance' => total_distance,
        'total_time' => total_time,
        'total_ascent' => total_ascent,
        'total_descent' => total_ascent,
        'average_heart_rate' => average_heart_rate,
        'maximum_heart_rate' => maximum_heart_rate,
        'average_power' => average_power,
        'total_calories' => total_calories
      }
    end

    def add_point(args)
      point = GeoPoint.new(args)
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

      if point.has_location? && @end_point.has_location?
        if point.speed.nil? || point.speed.nan?
          if point.time && @end_point.time
            point.speed = Maths.haversine_distance(@end_point, point) / (point.time - @end_point.time)
          end
        end
      end

      @end_point = point if point.has_location?
      get_points << point
    end

    def smooth_hr_lowess!(params)
      params[:delta] ||= 0

      valid = ->(p) { !p.heart_rate.nil? && (p.heart_rate.integer? || !p.heart_rate.nan?) }
      valid_points = get_points.select(&valid)
      return if valid_points.empty?

      lw_hr = Lowess.lowess(valid_points.map { |p| Lowess::Point.new(p.time, p.heart_rate) }, **params)
      i = 0
      get_points.map! do |p|
        if valid.call(p)
          p.heart_rate = lw_hr[i].y
          i += 1
        end
        p
      end
    end

    def smooth_speed_lowess!(params)
      params[:delta] ||= 0

      valid = ->(p) { !p.speed.nil? && !p.time.nil? && !p.speed.nan? }
      valid_points = get_points.select(&valid)
      return if valid_points.empty?

      lw_speed = Lowess.lowess(valid_points.map { |p| Lowess::Point.new(p.time, p.speed) }, **params)
      i = 0
      get_points.map! do |p|
        if valid.call(p)
          p.speed = lw_speed[i].y
          i += 1
        end
        p
      end
    end

    # See https://www.rdocumentation.org/packages/gplots/versions/3.0.1/topics/lowess
    def smooth!(speed_params, hr_params)
      smooth_speed_lowess! speed_params if speed_params != {}

      smooth_hr_lowess! hr_params if hr_params != {}
    end

    # A custom HR smoothing method for man-powered endurance sports.
    #
    # We smooth the heart rate with a LOWESS filter with parameters chosen to
    # match the maximum rate of rising heartrate on start of exercise and falling
    # heartrate in rest periods.
    # We there need a time scale which is used to compute linear fits, i.e. a
    # time where a linear fit is a good approximation.
    # Based on analysis of a set of HR traces this should be no more than 5 seconds.
    # We thus calculate the factor f for a lowess smoothing based on the total
    # time and a fraction of 10 seconds.
    def smooth_hr_endurance_sports!(ref_time = 10.0)
      smooth_hr_lowess! f: ref_time / total_time, iter: 4
    end

    # A simple speed smoother for on-water rowing taking into account typical acceleration times
    # Stroke rates are between 10 and 30 strokes per minute. At the very least
    # we need to average over several strokes to average out the front-and-aft
    # motion of wrist-worn gps receivers.
    # This means a minimum smoothing time of 15 seconds for 5 strokes at 20 SPM.
    #
    # (Ideally we would use the characteristics of a rowing
    # boat (such as turn radius etc) to do Kalman filtering of the GPS coordinates
    # and then calculate a new speed from this)
    def smooth_speed_rowing!(ref_time = 15.0)
      smooth_speed_endurance_sports! ref_time
    end

    # assume that for most endurance sports the speed is approximately linear
    # at 10-second timescales
    def smooth_speed_endurance_sports!(ref_time = 10.0)
      smooth_speed_lowess! f: ref_time / total_time, iter: 4
    end

    def add_lap(args)
      get_laps << Lap.new(args)
    end

    def process_elevation_delta(last_point, new_point)
      if last_point && new_point && last_point.elevation && new_point.elevation
        delta = new_point.elevation - last_point.elevation
        @_total_ascent = total_ascent + (delta > 0 ? delta : 0)
        @_total_descent = total_descent - (delta < 0 ? delta : 0)
      end
    end

    def started_at
      return @started_at if @started_at

      @start_point&.time
    end

    def ended_at
      @end_point&.time
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
      @total_distance&.round
    end

    # Public : Measure of how hilly the route is. Measured as total ascent (m) / distance (km)
    #
    # Returns Float measure
    def hilliness
      total_distance > 0 ? (total_ascent * 1000 / total_distance) : 0
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
      points.map(&:speed).compact.max || 0.0
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
      points.map(&:speed).compact.min || 0.0
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
      points.map(&:heart_rate).compact.max || 0
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
      points.map(&:heart_rate).compact.min || 0
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
      points.map(&:elevation).compact.max || 0
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
      points.map(&:elevation).compact.min || 0
    end

    # Public: Get total calories for whole GeoRoute.
    #
    # Examples
    #   @route.total_calories
    #   # => 10
    #
    # Returns Integer calories, 0 if no calories on laps or no laps on the route.
    def total_calories
      laps.map(&:calories).inject { |sum, l| sum + l } || 0
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
