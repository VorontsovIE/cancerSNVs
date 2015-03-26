require_relative 'sequence'

class SequenceWithSNP
  attr_reader :left, :allele_variants, :right, :name
  def initialize(left, allele_variants, right, options = {})
    # raise "SequenceWithSNP left part is invalid: #{left}" unless Sequence.valid_sequence?(left)
    # raise "SequenceWithSNP right part is invalid: #{right}" unless Sequence.valid_sequence?(right)
    # raise "SequenceWithSNP allele_variants are invalid: #{allele_variants}" unless allele_variants.map(&:to_s).all?{|letter| %w[A C G T N].include?(letter.upcase) }
    @left, @allele_variants, @right = left, allele_variants.map(&:to_s), right
    @name = options[:name] || (left + '_' + allele_variants.join('-') + '_' + right)
  end

  def self.from_string(sequence, options = {})
    left, mid, right = sequence.split(/[\[\]]/)
    allele_variants = mid.split('/')
    SequenceWithSNP.new(left, allele_variants, right, options)
  end

  def length
    left.length + 1 + right.length
  end

  def revcomp
    SequenceWithSNP.new(Sequence.revcomp(right),
                        allele_variants.map{|letter| Sequence.complement(letter) },
                        Sequence.revcomp(left), name: name)
  end

  def pyrimidine_context?
    ['C', 'T'].include?(allele_variants.first.upcase)
  end

  def in_pyrimidine_context
    pyrimidine_context? ? self : self.revcomp
  end

  def ==(other)
    other.is_a?(SequenceWithSNP) && @left == other.left && @allele_variants == other.allele_variants && @right == other.right
  end

  def eql?(other)
    other.class.equal?(SequenceWithSNP) && @left.eql?(other.left) && @allele_variants.eql?(other.allele_variants) && @right.eql?(other.right)
  end

  def hash
    [@left, @allele_variants, @right].hash
  end

  # main variant (or allele_variant_number variant) context
  def context(before: 1, after: 1, allele_variant_number: 0)
    left[-before..-1] + allele_variants[allele_variant_number] + right[0, after]
  end

  def subsequence(before:, after:)
    SequenceWithSNP.new(left[-before..-1], allele_variants, right[0, after], name: name)
  end

  def sequence_variant(allele_variant_number)
    "#{left}#{allele_variants[allele_variant_number]}#{right}"
  end

  def shuffle_string(str)
    str.each_char.to_a.shuffle.join
  end
  private :shuffle_string

  def shuffle_flanks(name: nil)
    shuffled_left = shuffle_string(left[0..-2]) + left[-1] # preserve 1-bp context
    shuffled_right = right[0] + shuffle_string(right[1..-1])
    SequenceWithSNP.new(shuffled_left, allele_variants, shuffled_right, name: name)
  end

  def to_s
    "#{name}\t#{left}[#{allele_variants.join('/')}]#{right}"
  end

  def inspect
    to_s
  end
end
