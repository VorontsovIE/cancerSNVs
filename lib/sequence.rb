class Sequence
  attr_reader :sequence, :name
  def initialize(sequence, options = {})
    raise "Wrong sequence `#{sequence}`"  unless Sequence.valid_sequence?(sequence)
    @sequence = sequence
    @name = options[:name] || sequence
  end

  def length
    sequence.length
  end

  def revcomp
    Sequence.new(Sequence.revcomp(sequence), name: name)
  end

  def self.complement(sequence)
    sequence.tr('acgtnACGTN'.freeze, 'tgcanTGCAN'.freeze)
  end

  def self.revcomp(sequence)
    complement(sequence).reverse
  end

  def self.valid_sequence?(sequence)
    sequence.match /\A[acgt]+\z/i
  end
end
