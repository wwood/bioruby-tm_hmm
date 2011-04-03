#!/usr/bin/env ruby

require 'bio'
require 'bio-tm_hmm'
require 'optparse'

# If being run directly instead of being require'd,
# output one transmembrane per line, and
# indicate that a particular protein has no transmembrane domain
if $0 == __FILE__
  options = {
    :filter_in => false,
    :filter_out => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = [
      'Usage: tm_hmm_wrapper.rb [-fg] [fasta_filename]',
      "\tfasta file can also be piped in on STDIN.",
      "\twithout arguments, a description of the transmembranes is printed out for each input sequence"
    ]
    opts.on('-f','--filter-in','Print those sequences that have a transmembrane domain') do
      options[:filter_in] = true
    end
    opts.on('-o','--filter-out','Print those sequences that do NOT have a transmembrane domain') do
      options[:filter_out] = true
    end
  end
  o.parse!

  runner = Bio::TMHMM::TmHmmWrapper.new

  Bio::FlatFile.auto(ARGF).each do |seq|
    result = runner.calculate(seq.seq)
    name = seq.definition

    # Default output - a description of the TMDs for each input aaseq
    if options[:filter_in] == false and options[:filter_out] == false
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

    elsif options[:filter_in]
    	if result.has_domain?
    		puts seq
    	end
    elsif options[:filter_out]
    	unless result.has_domain?
    		puts seq
    	end
    end
  end
end
