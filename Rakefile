# frozen_string_literal: true

require 'rake'
require 'rspec/core/rake_task'
require_relative 'lib/broutes'

task default: :spec
RSpec::Core::RakeTask.new

task :parse, :input_file, :format do |_t, args|
  Broutes.logger.level = Logger::DEBUG
  puts args
  file = File.open(args[:input_file])
  route = Broutes.from_file(file, args[:format].to_sym)
  file.close
  puts route.summary
  puts route.to_hash
end

task :plot, :input_file, :format do |_t, args|
  Broutes.logger.level = Logger::DEBUG
  puts args
  file = File.open(args[:input_file])
  route = Broutes.from_file(file, args[:format].to_sym)
  route.smooth_hr_endurance_sports!
  route.smooth_speed_rowing! 20.0
  puts route.points.to_a.length
  route.drop_to_target_rate! 0.20
  puts route.points.to_a.length
  file.close
  puts route.summary

  x = route.points.map { |p| p.time - route.started_at }
  y1 = route.points.map { |p| p.speed ? p.speed * 3.6 : nil } # km/h
  y2 = route.points.map(&:heart_rate)

  require 'tempfile'
  f = Tempfile.new ['gnuplot', '.tsv']
  begin
    f.write x.zip(y1, y2).map { |r| r.join "\t" }.join("\n")

    require 'numo/gnuplot'
    Numo.gnuplot do
      debug_on
      send 'set datafile separator "\t"'
      set title: args[:input_file]
      set :y2tics
      set xlabel: 'time [s]'
      set ylabel: 'speed [km/h]'
      send "set y2label 'HR (bpm)'"

      send "plot '#{f.path}' u 1:2 w l t 'speed' axis x1y1, \
                 '#{f.path}' u 1:3 w l t 'HR' axis x1y2"
      pause 50
    end
  ensure
    f.close
    f.unlink
  end
end

def parse_dir(dirname)
  Dir.foreach(dirname) do |item|
    next if (item == '.') || (item == '..')

    full_path = dirname + '/' + item
    if item.end_with?('.fit') ||
       item.end_with?('.tcx') ||
       item.end_with?('.gpx')
      puts item
      file = File.open(full_path)
      type = item.split('.')[-1].to_sym
      type = :gpx_track if type == :gpx
      route = Broutes.from_file(file, type)
      file.close
      puts route.summary
    elsif File.directory?(full_path)
      parse_dir(full_path)
    end
  end
end

# Recursively parse a directory of files
task :parse_dir, :dirname do |_t, args|
  Broutes.logger.level = Logger::ERROR
  puts args
  parse_dir(args[:dirname])
end
