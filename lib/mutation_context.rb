require 'set'

class MutationContext
  NUCLEOTIDE_SET = Set.new(%w[A C G T])
  attr_reader :before, :center, :after, :should_raise
  def initialize(before, center, after, should_raise: false)
    @before = before.upcase
    @center = center.upcase
    @after = after.upcase
    @should_raise = should_raise # raise in unknown letters in context on matching
    raise 'Context nucleotides should be one of: A,C,G,T or N ' unless %w[A C G T N].include?(@before) && %w[A C G T N].include?(@center) && %w[A C G T N].include?(@after)
    define_matcher!
  end

  def self.from_string(str, should_raise: false)
    raise 'Context should have length 3' unless str.length == 3
    self.new(str[0], str[1], str[2], should_raise: should_raise)
  end

  def revcomp
    MutationContext.new(Sequence.complement(after), Sequence.complement(center), Sequence.complement(before), should_raise: should_raise)
  end

  private def define_matcher!
    check_before = (@before == 'N') ? nil : "('#{@before}' == str[0])"
    check_center = (@center == 'N') ? nil : "('#{@center}' == str[1])"
    check_after  = (@after == 'N')  ? nil : "('#{@after}'  == str[2])"
    check = [check_before, check_center, check_after].compact
    result = check.empty? ? 'true' : check.join(' && ')
    if should_raise
      raising = <<-EOS
        raise "Undefined letter in string `\#{str}`"  unless NUCLEOTIDE_SET.include?(str[0]) && NUCLEOTIDE_SET.include?(str[1]) && NUCLEOTIDE_SET.include?(str[2])
      EOS
    else
      raising = ''
    end
    meth_body = <<-EOS
      def match?(str)
        #{raising}
        #{result}
      end
    EOS
    instance_eval(meth_body)
  end

  def to_s; "MutationContext<#{before}#{center}#{after}>"; end
  def inspect; "MutationContext<#{before}#{center}#{after}>"; end

  def hash; [before, center, after].hash; end;
  def eql?(other); other.class.equal?(self.class)  &&  before == other.before  &&  center == other.center  &&  after == other.after; end
  def ==(other); eql?(other); end    
end
