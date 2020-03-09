# frozen_string_literal: true

module Broutes
  require 'logger'
  require_relative 'broutes/geo_route'
  require_relative 'broutes/geo_point'
  require_relative 'broutes/lap'
  require_relative 'broutes/maths'
  require_relative 'broutes/formats'

  RAD_PER_DEG = 0.017453293 #  PI/180
  EARTH_RADIUS = 6_371_000 # m

  def self.from_file(file, format, speed_params = {}, hr_params = {})
    raise "unable to interpret format #{format}" unless processor = Formats::Factory.new.get(format)

    Broutes.logger.debug { "found processor #{processor} for #{file}" }
    route = GeoRoute.new
    processor.load(file, route)
    route.smooth!(speed_params, hr_params)
    route
  end

  def self.from_hash(h)
    GeoRoute.new h
  end

  class << self
    attr_writer :logger

    def logger
      @logger ||= Logger.new(STDOUT)
    end
  end
end
