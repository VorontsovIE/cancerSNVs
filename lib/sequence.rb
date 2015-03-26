class Sequence
  attr_reader :sequence, :name
  def initialize(sequence, options = {})
    # raise "Wrong sequence `#{sequence}`"  unless Sequence.valid_sequence?(sequence)
    @sequence = sequence
    @name = options[:name] || sequence
  end

  def length
    sequence.length
  end

  def revcomp
    Sequence.new(Sequence.revcomp(sequence), name: name)
  end

  def ==(other)
    other.is_a?(Sequence) && @sequence == other.sequence
  end

  def eql?(other)
    other.class.equal?(Sequence) && @sequence.eql?(other.sequence)
  end

  def hash
    @sequence.hash
  end

  def self.complement(sequence)
    sequence.tr('acgtnACGTN'.freeze, 'tgcanTGCAN'.freeze)
  end

  def self.revcomp(sequence)
    complement(sequence).reverse
  end

  def self.valid_sequence?(sequence)
    sequence.match /\A[acgtn]+\z/i
  end
end
