module GenomeReader
  # Genome loaded into memory
  class MemoryReader
    # chromosome names are Symbols, not Strings
    def self.load_from_disk(genome_folder,
                            chromosome_name_matcher: /^(?<chromosome>\w+)\.plain$/)
      chromosome_sequences = Dir.glob(File.join(genome_folder, '*')).map{|filename|
        chromosome_name = File.basename(filename)[ chromosome_name_matcher, :chromosome ]
        [chromosome_name.to_sym, File.read(filename)]
      }.to_h

      self.new(chromosome_sequences)
    end

    def initialize(chromosome_sequences)
      @chromosome_sequences = chromosome_sequences
    end

    # make sure, chromosome is a Symbol
    def chromosome_exist?(chromosome)
      @chromosome_sequences.has_key?(chromosome)
    end

    def chromosome_names
      @chromosome_sequences.keys
    end

    def read_sequence(chromosome, coordinate_system, from, to)
      pos = coordinate_system.seek_position(from)
      len = coordinate_system.length(from, to)

      if len.is_a?(Integer) || len.to_f.finite? # finite length
        chromosome[pos, len]
      else
        chromosome[pos..-1]
      end
    end

    def to_s
      chromosomes_text = chromosome_names.sort.first(5).join(', ') + (chromosome_names.size > 5 ? ', ...' : '')
      "GenomeReader::MemoryReader<#{chromosome_sequences.size} chromosomes: #{chromosomes_text}>"
    end

    def inspect
      to_s
    end
  end
end
