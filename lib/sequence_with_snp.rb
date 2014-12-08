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
                        Sequence.revcomp(left))
  end

  # main variant (or allele_variant_number variant) context
  def context(before: 1, after: 1, allele_variant_number: 0)
    left[-before..-1] + allele_variants[allele_variant_number] + right[0, after]
  end
end
