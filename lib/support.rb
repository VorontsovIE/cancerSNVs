def count_each_element(values)
  result = {}
  values.each do |value|
    if result.has_key?(value)
        result[value] += 1
    else
        result[value] = 1
    end
  end
  result
end

module BreastCancerSNPs
  class OutputConfigurator
    attr_reader :filename_result, :motif_names
    def initialize(filename_result, motif_names)
      @filename_result, @motif_names = filename_result, motif_names
      Dir.mkdir dirname unless Dir.exist?(dirname)
    end

    def dirname
      File.dirname(filename_result)
    end

    def basename
      File.basename(filename_result)
    end

    def filename_w_prefix(prefix)
      File.join(dirname, prefix + basename)
    end

    def output_counts_for_each_motif(output_filename, motif_counts)
      File.open(output_filename, 'w') do |fw|
        motif_names.each do |motif_name|
          fw.puts "#{motif_name}\t#{motif_counts.fetch(motif_name, 0)}"
        end
      end
    end

    def motif_statistics(file_prefix, motif_counts)
      output_counts_for_each_motif(filename_w_prefix(file_prefix), motif_counts)
    end
  end
end
