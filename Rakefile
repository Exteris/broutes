require 'rake'
require 'rspec/core/rake_task'
require_relative 'lib/broutes'

task :default => :spec
RSpec::Core::RakeTask.new

task :parse, :input_file, :format do |t, args|
  Broutes.logger.level = Logger::ERROR
  puts args
  file = File.open(args[:input_file])
  route = Broutes.from_file(file, args[:format].to_sym)
  file.close
  puts route.summary
end

def parse_dir(dirname)
  Dir.foreach(dirname) do |item|
    next if item == '.' or item == '..'
    full_path = dirname+'/'+item
    if (item.end_with?('.fit') or
        item.end_with?('.tcx') or
        item.end_with?('.gpx'))
      file = File.open(full_path)
      type = item.split('.')[-1].to_sym
      type = :gpx_track if type == :gpx
      route = Broutes.from_file(file, type)
      file.close
      puts item, route.summary
    elsif File.directory?(full_path)
      parse_dir(full_path)
    end
  end
end

# Recursively parse a directory of files
task :parse_dir, :dirname do |t, args|
  Broutes.logger.level = Logger::ERROR
  puts args
  parse_dir(args[:dirname])
end
