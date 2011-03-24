#!/usr/bin/env ruby

$:.unshift File.join(File.dirname(__FILE__),'..','vendor/plugins/ben_bioinformatics/lib')

require 'rubygems'
require 'transmembrane'
include Transmembrane
require 'tempfile'
gem 'rio'
require 'rio'

class TmHmmWrapper
  
  # Given an amino acid sequence, return a TransmembraneProtein
  # made up of the predicted transmembrane domains
  def calculate(sequence)
    rio(:tempdir) do |d|
      FileUtils.cd(d.to_s) do
        Tempfile.open('tmhmmin') { |tempfilein|
          # Write a fasta to the tempfile
          tempfilein.puts '>wrapperSeq'
          tempfilein.puts "#{sequence}"
          tempfilein.close #required. Maybe because it doesn't flush otherwise?
      
          Tempfile.open('signalpout') {|out|
            result = system("tmhmm -short #{tempfilein.path} >#{out.path}")
        
            if !result
              raise Exception, "Running TMHMM program failed. See $? for details."
            end
        
        
            line = rio(out.path).readline
            return TmHmmResult.create_from_short_line(line)
          }
        }
      end
    end
  end
end


class TmHmmResult
  attr_reader :domains
  
  # initialise with the output line of a 
  # eg. 
  #PFF0290w	len=293	ExpAA=145.77	First60=20.51	PredHel=7	Topology=o39-61i101-120o140-162i169-186o196-218i230-252o262-284i
  def self.create_from_short_line(line)
    protein = OrientedTransmembraneDomainProtein.new
    
    splits = line.strip.split("\t")
    if splits.length != 6
      raise Exception, "Incorrectly parsed short line from TMHMM: #{line}"
    end
    
    substrate = splits[5]
    if substrate.gsub!(/^Topology\=[io]/,'').nil?
      raise Exception, "Badly parsed Topology hit: #{substrate}"
    end
    
    matches = substrate.match('^(\d+?)\-')
    if !matches
      return protein #no transmembrane domains predicted
    end
    
    # eat the string from the beginning adding the transmembrane domains
    prev = matches[1]
    substrate.gsub!(/^(\d+?)-/,'')
    # match all the middle bits
    reg = /^(\d+?)([io])(\d+?)\-/
    while matches =substrate.match(reg)
      tmd = OrientedTransmembraneDomain.new
      tmd.start = prev.to_i
      tmd.stop = matches[1].to_i
      tmd.orientation = parse_orientation_from_last_location(matches[2])
      protein.push tmd
      
      prev = matches[3]
      substrate.gsub!(reg, '')
    end
    #match the last bit
    if !(matches = substrate.match('(\d+?)([io])$'))
      raise Exception, "Failed to parse the last bit of: #{substrate}"
    end
    tmd = OrientedTransmembraneDomain.new
    tmd.start = prev.to_i
    tmd.stop = matches[1].to_i
    tmd.orientation = parse_orientation_from_last_location(matches[2])
    protein.push tmd
    
    return protein
  end
  
  def self.parse_orientation_from_last_location(last_location)
    case last_location
    when 'i'  
      return OrientedTransmembraneDomain::OUTSIDE_IN
    when 'o'
      return OrientedTransmembraneDomain::INSIDE_OUT
    else
      raise Exception, "Badly parsed topology hit due to orientation character: #{substrate}"
    end
  end
end


# If being run directly instead of being require'd, 
# output one transmembrane per line, and
# indicate that a particular protein has no transmembrane domain
if $0 == __FILE__
  require 'bio'
  
  runner = TmHmmWrapper.new
  
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
