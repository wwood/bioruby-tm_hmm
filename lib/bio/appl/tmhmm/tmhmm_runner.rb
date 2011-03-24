require 'tempfile'

module Bio
  class TMHMM
    class TmHmmWrapper
      # Given an amino acid sequence, return a TransmembraneProtein
      # made up of the predicted transmembrane domains
      def calculate(sequence)
        Bio::Command.mktmpdir do |d|
          line = nil
          Bio::Command.call_command(['tmhmm','-short'], :chdir => d) do |io|
            io.puts '>wrapperSeq'
            io.puts sequence
            io.close_write
            line = io.readline
          end
          
          if line.nil?
            raise Exception, "Error running locally installed TMHMM program 'tmhmm'. Is it properly installed?"
          end
          
          return TmHmmResult.create_from_short_line(line)
        end
      end
    end
    
    class TmHmmResult
      attr_reader :domains
      
      # initialise with the output line of a 
      # eg. 
      #PFF0290w len=293 ExpAA=145.77  First60=20.51 PredHel=7 Topology=o39-61i101-120o140-162i169-186o196-218i230-252o262-284i
      def self.create_from_short_line(line)
        protein = Bio::Transmembrane::OrientedTransmembraneDomainProtein.new
        
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
          tmd = Bio::Transmembrane::OrientedTransmembraneDomain.new
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
        tmd = Bio::Transmembrane::OrientedTransmembraneDomain.new
        tmd.start = prev.to_i
        tmd.stop = matches[1].to_i
        tmd.orientation = parse_orientation_from_last_location(matches[2])
        protein.push tmd
        
        return protein
      end
      
      def self.parse_orientation_from_last_location(last_location)
        case last_location
          when 'i'  
          return Bio::Transmembrane::OrientedTransmembraneDomain::OUTSIDE_IN
          when 'o'
          return Bio::Transmembrane::OrientedTransmembraneDomain::INSIDE_OUT
        else
          raise Exception, "Badly parsed topology hit due to orientation character: #{substrate}"
        end
      end
    end
  end
end