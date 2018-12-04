require 'active_support'
require 'fit_parser'

module Broutes::Formats
  class FitFile

    def load(file, route)
      fit_file = FitParser.load_file(file)
      Broutes.logger.info {"Started fit processing"}
      i = 0
      fit_file.records.select {|r| r.content && r.content.record_type == :record }.each do |r|
        begin
          pr = r.content
          data = { time: record_time(r) }
          data[:lat] = convert_position(pr.raw_position_lat) if pr.respond_to?(:raw_position_lat)
          data[:lon] = convert_position(pr.raw_position_long) if pr.respond_to?(:raw_position_long)
          data[:elevation] = pr.raw_altitude if pr.respond_to?(:raw_altitude) and pr.raw_altitude.to_i != 65535
          [:distance, :heart_rate, :power, :speed, :cadence, :temperature].each do |m| 
            data[m] = pr.send("raw_#{m}").to_f if pr.respond_to?("raw_#{m}") and not [255, 65535].include?(pr.send("raw_#{m}").to_i)
          end

          route.add_point(data)
          i += 1
        rescue => e
          Broutes.logger.debug {"#{e.message} for #{r}"}
        end
      end
      Broutes.logger.info {"Loaded #{i} data points"}
    end

    def convert_position(value)
      val = (8.381903171539307e-08 * value).round(5)
      if val == 180.0 # special value for some watches
        return nil
      else
        return val
      end
    end

    def record_time(record)
      utc_seconds = record.content.raw_timestamp
      utc_seconds += record.header.time_offset if record.header.compressed_timestamp? 
      Time.new(1989, 12, 31, 0, 0, 0, "+00:00") + utc_seconds #seconds since UTC 00:00 Dec 31 1989
    end

  end
end
