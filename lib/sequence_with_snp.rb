require_relative 'sequence'

class SequenceWithSNP
  attr_reader :left, :allele_variants, :right, :name
  def initialize(left, allele_variants, right, options = {})
    raise unless Sequence.valid_sequence?(left)
    raise unless Sequence.valid_sequence?(right)
    raise unless allele_variants.all?{|letter| %w[A C G T].include?(letter.upcase) }
    @left, @allele_variants, @right = left, allele_variants, right
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

  # main variant (or allele_variant_number variant) context
  def context(before: 1, after: 1, allele_variant_number: 0)
    left[-before..-1] + allele_variants[allele_variant_number] + right[0, after]
  end

  def shuffle_string(str)
    str.each_char.to_a.shuffle.join
  end
  private :shuffle_string

  def shuffle_flanks(name: nil)
    SequenceWithSNP.new(shuffle_string(left), allele_variants, shuffle_string(right), name: name)
  end

  def to_s
    "#{name}\t#{left}[#{allele_variants.join('/')}]#{right}"
  end
end
