# Attention! This file relies on experiment configuration!
# (because it loads motifs from collection folder which is set by relative path)
require_relative '../../experiment_configuration'
require 'set'
require 'bioinform'

module PerfectosAPE
    
  ResultShort = Struct.new( :line,
                            :variant_id, :motif_name,
                            :pvalue_1, :pvalue_2,
                            :pos_1, :orientation_1,
                            :pos_2, :orientation_2) do
    
    # it's impossible to get length from short site info, so we preload motifs
    class << self
      attr_reader :motif_lengths
    end
    @motif_lengths = {}
    Dir.glob(File.join(MOTIF_COLLECTION_FOLDER, '*.pwm')) do |motif_filename|
      motif = Bioinform::MotifModel::PWM.from_file(motif_filename)
      @motif_lengths[motif.name.to_sym] ||= motif.length
    end

    # "27610826_3 MAZ_f1   1.12e-04  9.60e-04  -7  +  -7  +"
    def self.from_string(line)
      variant_id, motif_name, pvalue_1, pvalue_2,  pos_1, orientation_1,  pos_2, orientation_2 = line.split("\t")
      self.new( line,
                variant_id, motif_name.to_sym,
                pvalue_1.to_f, pvalue_2.to_f,
                pos_1.to_i,  (orientation_1 == '+') ? :direct : :revcomp,
                pos_2.to_i,  (orientation_2 == '+') ? :direct : :revcomp )
    end

    def normalized_snv_name
      variant_id.split("_", 2)[0]
    end

    def length
      PerfectosAPE::ResultShort.motif_lengths[motif_name]
    end

    def seq_1_five_flank_length
      -pos_1
    end

    def seq_1_three_flank_length
      length - 1 + pos_1
    end

    def fold_change
      pvalue_1 / pvalue_2
    end

    def site_before_substitution?(pvalue_cutoff: 0.0005)
      pvalue_1 <= pvalue_cutoff
    end

    def site_after_substitution?(pvalue_cutoff: 0.0005)
      pvalue_2 <= pvalue_cutoff
    end

    def disrupted?(fold_change_cutoff: 5)
      fold_change <= (1.0 / fold_change_cutoff)
    end

    def emerged?(fold_change_cutoff: 5)
      fold_change >= fold_change_cutoff
    end

    # glues corresponding motif site positions on direct and reverse strands
    def snv_position_in_site_1_pwm
      if orientation_1 == :direct
        - pos_1
      else
        pos_1 + length - 1
      end
    end

    def substitution_in_core?
      pos = snv_position_in_site_1_pwm
      pos >= 0 && pos < length
    end

    def substitution_in_flank?
      ! substitution_in_core?
    end

    def self.each_in_stream(stream, &block)
      stream.each_line.lazy.reject{|line|
        line.start_with?('#')
      }.map{|line|
        self.from_string(line)
      }.each(&block)
    end

    def self.each_in_file(all_mutations_filename, &block)
      return enum_for(:each_in_file, all_mutations_filename).lazy  unless block_given?
      File.open(all_mutations_filename) do |f|
        each_in_stream(f, &block)
      end
    end
  end
end
