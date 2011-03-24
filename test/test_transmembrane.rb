# 
# To change this template, choose Tools | Templates
# and open the template in the editor.
 

$:.unshift File.join(File.dirname(__FILE__),'..','lib')

require 'test/unit'
require 'transmembrane'

module Transmembrane
  class TransmembraneTest < Test::Unit::TestCase
    def test_confidenced_transmembrane_domain
      one = ConfidencedTransmembraneDomain.new
      two = ConfidencedTransmembraneDomain.new
      assert_equal one, two
    end
    
    def test_sequence_offsets
      aaseq = 'AAAAAANG' #8 aa long
      d = TransmembraneDomainDefinition.new
      d.start = 6
      d.stop = 8
      assert_equal 'ANG', d.sequence(aaseq)
      
      assert_equal 'AANG', d.sequence('AAAAAANG', -1, 0)
      assert_equal 'AANG', d.sequence('AAAAAANG', -1, 1) #overhang
      assert_equal 'AAN', d.sequence('AAAAAANG', -1, -1) #overhang over the cterm
      
      d.start = 1
      d.stop = 5
      assert_equal 'AAAAA', d.sequence('AAAAAANG', -2, 0) #overhang over the nterm
      assert_equal 'AAAAAANG', d.sequence('AAAAAANG', -2, 15) #overhang over the nterm and cterm
    end
  end

  class TransmembraneProteinTest < Test::Unit::TestCase
    def test_simple
      protein = TransmembraneProtein.new
      tmd = TransmembraneDomainDefinition.new
      tmd.start = 8
      tmd.stop = 9
      protein.push tmd

      tmd = TransmembraneDomainDefinition.new
      tmd.start = 8
      tmd.stop = 10
      protein.push tmd

      assert_equal 2, protein.minimum_length
      assert_equal 2.5, protein.average_length
    end

    def test_empty
      protein = TransmembraneProtein.new
      assert protein.transmembrane_domains.empty?
      assert_equal false, protein.has_domain?
    end

    def test_overlaps
      p1 = TransmembraneProtein.new
      p2 = TransmembraneProtein.new
      tmd1 = TransmembraneDomainDefinition.new(5,10)
      tmd2 = TransmembraneDomainDefinition.new(5,6)
      p1.transmembrane_domains = [tmd1]
      p2.transmembrane_domains = [tmd2]
      assert_equal [[tmd1, tmd2]], p1.overlaps(p2)

      p2.transmembrane_domains = [tmd1,tmd2]
      assert_equal [[tmd1, tmd1],[tmd1, tmd2]], p1.overlaps(p2)

      tmd3 = TransmembraneDomainDefinition.new(500,600)
      p2.transmembrane_domains = [tmd3]
      assert_equal [], p1.overlaps(p2)
    end

    def test_best_overlap
      p1 = TransmembraneProtein.new
      p2 = TransmembraneProtein.new
      tmd1 = TransmembraneDomainDefinition.new(5,10)
      tmd2 = TransmembraneDomainDefinition.new(5,6)
      tmd3 = TransmembraneDomainDefinition.new(11,22)
      p1.transmembrane_domains = [tmd1]
      p2.transmembrane_domains = [tmd2]
      assert_equal [tmd1, tmd2], p1.best_overlap(p2)

      p2.transmembrane_domains = [tmd1,tmd2]
      assert_equal [tmd1, tmd1], p1.best_overlap(p2)

      p2.transmembrane_domains = [tmd3]
      assert_equal nil, p1.best_overlap(p2)
    end

    def test_each
      expected = [10,6,22]
      p = TransmembraneProtein.new
      tmd1 = TransmembraneDomainDefinition.new(5,10)
      tmd2 = TransmembraneDomainDefinition.new(5,6)
      tmd3 = TransmembraneDomainDefinition.new(11,22)
      p.transmembrane_domains = [tmd1,tmd2,tmd3]
      p.each_with_index do |tmd, i|
        assert_equal expected[i], tmd.stop
      end
    end
  end

  class TransmembraneProteinTest < Test::Unit::TestCase
    def test_overlap_tmd
      tmd1 = TransmembraneDomainDefinition.new(3,6)
      tmd2 = TransmembraneDomainDefinition.new(3,6)
      assert_equal((3..6), tmd1.intersection(tmd2))
      assert_equal 4, tmd1.overlap_length(tmd2)

      tmd2 = TransmembraneDomainDefinition.new(7,8)
      assert_equal 0, tmd1.overlap_length(tmd2)

      tmd2 = TransmembraneDomainDefinition.new(5,8)
      assert_equal 2, tmd1.overlap_length(tmd2)

      tmd2 = TransmembraneDomainDefinition.new(1,3)
      assert_equal 1, tmd1.overlap_length(tmd2)
      assert_equal((3..3), tmd1.intersection(tmd2))
    end
  end
end
