# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-
# stub: bio-tm_hmm 0.2.4 ruby lib

Gem::Specification.new do |s|
  s.name = "bio-tm_hmm".freeze
  s.version = "0.2.4"

  s.required_rubygems_version = Gem::Requirement.new(">= 0".freeze) if s.respond_to? :required_rubygems_version=
  s.require_paths = ["lib".freeze]
  s.authors = ["Ben J. Woodcroft".freeze]
  s.date = "2017-04-06"
  s.description = "A bioruby plugin for interaction with the transmembrane predictor TMHMM".freeze
  s.email = "donttrustben@gmail.com".freeze
  s.executables = ["bio-tm_hmm".freeze]
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.md"
  ]
  s.files = [
    ".document",
    ".travis.yml",
    "Gemfile",
    "LICENSE.txt",
    "README.md",
    "Rakefile",
    "VERSION",
    "bin/bio-tm_hmm",
    "bio-tm_hmm.gemspec",
    "lib/bio-tm_hmm.rb",
    "lib/bio/appl/tmhmm/tmhmm_runner.rb",
    "lib/bio/transmembrane.rb",
    "test/data/falciparum1.fa",
    "test/helper.rb",
    "test/test_tm_hmm_parser.rb",
    "test/test_tm_hmm_wrapper.rb",
    "test/test_transmembrane.rb"
  ]
  s.homepage = "http://github.com/wwood/bioruby-tm_hmm".freeze
  s.licenses = ["MIT".freeze]
  s.rubygems_version = "2.5.2".freeze
  s.summary = "A bioruby plugin for interaction with the transmembrane predictor TMHMM".freeze

  if s.respond_to? :specification_version then
    s.specification_version = 4

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<bio>.freeze, [">= 1.4.2"])
      s.add_development_dependency(%q<test-unit>.freeze, ["~> 3.2"])
      s.add_development_dependency(%q<shoulda>.freeze, [">= 0"])
      s.add_development_dependency(%q<rdoc>.freeze, ["~> 3.12"])
      s.add_development_dependency(%q<jeweler>.freeze, ["~> 2.3"])
      s.add_development_dependency(%q<bundler>.freeze, [">= 1.0.21"])
    else
      s.add_dependency(%q<bio>.freeze, [">= 1.4.2"])
      s.add_dependency(%q<test-unit>.freeze, ["~> 3.2"])
      s.add_dependency(%q<shoulda>.freeze, [">= 0"])
      s.add_dependency(%q<rdoc>.freeze, ["~> 3.12"])
      s.add_dependency(%q<jeweler>.freeze, ["~> 2.3"])
      s.add_dependency(%q<bundler>.freeze, [">= 1.0.21"])
    end
  else
    s.add_dependency(%q<bio>.freeze, [">= 1.4.2"])
    s.add_dependency(%q<test-unit>.freeze, ["~> 3.2"])
    s.add_dependency(%q<shoulda>.freeze, [">= 0"])
    s.add_dependency(%q<rdoc>.freeze, ["~> 3.12"])
    s.add_dependency(%q<jeweler>.freeze, ["~> 2.3"])
    s.add_dependency(%q<bundler>.freeze, [">= 1.0.21"])
  end
end

