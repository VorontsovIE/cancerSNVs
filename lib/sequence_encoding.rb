# Encoding sequence --> numeric array is a projection into indices but it shouldn't be reversible
# (sequence isn't necessary fully recoverable (i.e. letter case can be lost)
# Indices are now supposed to start with zero due to +alphabet_length+ definition
class SequenceEncoder
  # alphabet_length = 5
  # letter_to_index = {'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 4, 'a' => 0, 'c' => 1, 'g' => 2, 't' => 3, 'n' => 4}
  # index_to_letter = {0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', 4 => 'N'}

  attr_reader :letter_to_index, :index_to_letter, :alphabet_length

  # letter_to_index: {'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 4,
  #                   'a' => 0, 'c' => 1, 'g' => 2, 't' => 3, 'n' => 4}
  # index_to_letter: {0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', 4 => 'N'}
  def initialize(letter_to_index:, index_to_letter:)
    @letter_to_index = letter_to_index
    @index_to_letter = index_to_letter
    @alphabet_length = index_to_letter.size
  end

  # convert sequence to numerical code
  def encode_sequence(seq)
    Array.new(seq.length){|i| @letter_to_index[seq[i]] }
  end

  # Converts numeral index to sequence
  # +code_length+ is necessary while numerical code can contain zeros because it 
  # makes it impossible to calculate sequence length
  def decode_sequence(ind, code_length:)
    letters = []
    code_length.times do
      letters.unshift(@index_to_letter[ind % @alphabet_length])
      ind /= @alphabet_length
    end
    letters.join
  end

  def self.default_encoder
    SequenceEncoder.new(letter_to_index: {'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 4,
                                          'a' => 0, 'c' => 1, 'g' => 2, 't' => 3, 'n' => 4},
                        index_to_letter: {0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', 4 => 'N'})
  end

end
