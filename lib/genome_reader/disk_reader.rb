module GenomeReader
  # Each call to a genome loads necessary part from disk
  class DiskReader
    def initialize(genome_folder,
                  chromosome_file_by_name: ->(chr){ "#{chr}.plain" },
                  chromosome_name_matcher: /^(?<chromosome>\w+)\.plain$/)
      @genome_folder = genome_folder
      @chromosome_file_by_name = chromosome_file_by_name
      @chromosome_name_matcher = chromosome_name_matcher
    end

    def chromosome_filename(chromosome)
      File.join(@genome_folder, @chromosome_file_by_name.call(chromosome))
    end
    private :chromosome_filename

    def chromosome_exist?(chromosome)
      File.exist?(chromosome_filename.exist?(chromosome))
    end

    # chromosome names are Symbols, not Strings
    def chromosome_names
      Dir.glob(File.join(@genome_folder, '*')).map{|filename|
        File.basename(filename)[ @chromosome_name_matcher, :chromosome ]
      }.compact.map(&:to_sym)
    end

    def read_sequence(chromosome, coordinate_system, from, to)
      pos = coordinate_system.seek_position(from)
      len = coordinate_system.length(from, to)
      File.open(chromosome_filename(chromosome)){|f|
        f.seek(pos)
        if len.is_a?(Integer) || len.to_f.finite? # finite length
          seq = f.read(len)
        else
          seq = f.read
        end
        seq.b # ASCII is much faster than Unicode
      }
    end

    def read_chromosome(chromosome)
      File.read(chromosome_filename(chromosome)).b # ASCII is much faster than Unicode
    end

    def to_s
      chromosomes_text = chromosome_names.sort.first(5).join(', ') + (chromosome_names.size > 5 ? ', ...' : '')
      "GenomeReader::MemoryReader<#{chromosome_names.size} chromosomes: #{chromosomes_text}>"
    end

    def inspect
      to_s
    end
  end
end
