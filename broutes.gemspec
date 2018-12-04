$:.push File.expand_path("../lib", __FILE__)
require "broutes/version"

Gem::Specification.new do |s|
  s.name        = "broutes"
  s.version     = Broutes::VERSION.dup
  s.platform    = Gem::Platform::RUBY
  s.summary     = "Utilities for parsing and creating geo route files"
  s.email       = "adam.bird@gmail.com"
  s.homepage    = "http://github.com/adambird/broutes"
  s.description = "Utilities for parsing and creating geo route files"
  s.authors     = ['Adam Bird']

  s.files         = Dir["lib/**/*"]
  s.test_files    = Dir["spec/**/*"]
  s.require_paths = ["lib"]

  s.add_dependency('nokogiri', '~> 1.5.9')
  s.add_dependency('garmin-fit', '~> 0.0.2')
  s.add_dependency('bindata', '~> 1.4.5')
  s.add_dependency('lowess', '~> 1.0')
end
