require 'spec_helper'

describe Broutes do
  describe ".from_file" do
    before(:all) do
      @file = open_file('single_lap_gpx_track.gpx')
      @route = Broutes.from_file(@file, 'single_lap_gpx_track.gpx')
    end

    it "sets the start point lat" do
      @route.start_point.lat.should eq(52.9552055)
    end
    it "sets the start point lon" do
      @route.start_point.lon.should eq(-1.1558583)
    end
    it "sets the total distance" do
      @route.total_distance.should eq(7088)
    end
    it "sets the total ascent" do
      @route.total_ascent.round.should eq(34)
    end
    it "sets the total descent" do
      @route.total_descent.round.should eq(37)
    end
  end

  describe ".from_file_recalc_speed" do
  end

  describe ".from_fit_file" do
    before(:all) do
      @file = open_file('sample_two.fit')
      @route = Broutes.from_file(@file, 'sample_two.fit')
    end
    it "sets the total distance" do
      @route.total_distance.should eq(10644)
    end
    it "has reasonable bounds for the speed" do
      expect(@route.points.map(&:speed).compact.max).to be < 30
    end
    it "sets the speed correctly" do
      # this doesn't handle nonuniform time correctly
      speeds = @route.points.map(&:speed).compact
      mean_speed = speeds.sum / speeds.length

      ref_speed = @route.total_distance / @route.total_time
      expect(mean_speed).to be_between(ref_speed*0.85, ref_speed*1.15)
    end
  end
end
