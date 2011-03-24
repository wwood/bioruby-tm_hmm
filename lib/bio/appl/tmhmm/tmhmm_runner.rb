module Bio
  class TMHMM
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