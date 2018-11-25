require 'spec_helper'
require 'benchmark'

describe Formats::FitFile do
  describe "#load" do
    it "loads within 1 seconds" do

      time = Benchmark.measure do
        5.times do
          @file = open_file('long_track.fit')
          @route = GeoRoute.new
          @target = Formats::FitFile.new
          @target.load(@file, @route)
        end
      end
      

      time.total.should be <= 1
    end
  end
end
