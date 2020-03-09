# frozen_string_literal: true

require 'time'

module Broutes
  class GeoPoint
    attr_accessor :lat, :lon, :elevation, :distance, :heart_rate, :power, :speed, :cadence, :temperature
    attr_reader :time

    def initialize(args = {})
      args.each_pair { |key, value| send("#{key}=", value) if respond_to?("#{key}=") }
    end

    def self.from_hash(h)
      GeoPoint.new(h)
    end

    def has_location?
      lat && lon
    end

    def has_hr?
      !heart_rate.nil? && (heart_rate.integer? || !heart_rate.nan?)
    end

    def has_speed?
      !speed.nil? && !speed.nan?
    end

    def time=(value)
      @time = if value.is_a?(String)
                DateTime.parse(value).to_time
              else
                value
              end
    end

    def ==(other)
      lat == other.lat &&
        lon == other.lon &&
        elevation == other.elevation &&
        distance == other.distance &&
        time == other.time &&
        heart_rate == other.heart_rate &&
        power == other.power &&
        speed == other.speed &&
        cadence == other.cadence &&
        temperature == other.temperature
    end

    def to_hash
      h = {}
      h['lat'] = lat if lat
      h['lon'] = lon if lon
      h['elevation'] = elevation if elevation
      h['distance'] = distance if distance
      h['time'] = time if time
      h['heart_rate'] = heart_rate if heart_rate
      h['power'] = power if power
      h['speed'] = speed if speed
      h['cadence'] = cadence if cadence
      h['temperature'] = temperature if temperature
      h
    end
  end
end
