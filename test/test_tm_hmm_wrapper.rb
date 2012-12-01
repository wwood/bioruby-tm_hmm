require 'bio'
require 'helper'

class TmHmmWrapperTest < Test::Unit::TestCase
  include Bio::Transmembrane
  def test_wrapper
    prog = Bio::TMHMM::TmHmmWrapper.new
    seq = Bio::FlatFile.auto(File.join(File.dirname(__FILE__),'data','falciparum1.fa')).next_entry
    tmp = prog.calculate(seq.seq)
    assert tmp
    assert_equal false, tmp.has_domain?
  end
end
