#!/usr/bin/env ruby

require 'bio'
require 'bio-tm_hmm'

# If being run directly instead of being require'd, 
# output one transmembrane per line, and
# indicate that a particular protein has no transmembrane domain
if $0 == __FILE__
  require 'bio'
  
  runner = Bio::TMHMM::TmHmmWrapper.new
  
  Bio::FlatFile.auto(ARGF).each do |seq|
    result = runner.calculate(seq.seq)
    name = seq.definition
    
    if result.has_domain?
      # At least one TMD found. Output each on a separate line
      result.transmembrane_domains.each do |tmd|
        puts [
          name,
          result.transmembrane_type,
          tmd.start,
          tmd.stop,
          tmd.orientation
        ].join("\t")
      end
    else
      puts [
        name,
        'No Transmembrane Domain Found'
      ].join("\t")
    end
  end
end
