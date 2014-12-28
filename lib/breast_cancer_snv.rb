require 'set'
require 'interval_notation'
require_relative 'support'

BreastCancerSNV = Struct.new( :variant_id,
                              :sample_id, :chr, :position, :genome_build,

                              :ref_base_plus_strand, :mutant_base_plus_strand,

                              :five_prime_flanking_sequence_in_pyrimidine_context,
                              :ref_base_pyrimidine_context,
                              :mutant_base_pyrimidine_context,
                              :three_prime_flanking_sequence_in_pyrimidine_context,

                              :strand_of_mutation_in_pyrimidine_context,
                              :gene, :gene_id, :ccds_id, :transcript_id,

                              :gene_type, :mut_types,

                              :mRNA_mut_syntax, :cds_mut_syntax, :aa_mut_syntax,
                              :current_conf_status, :validation_platform ) do

  # '27607140	PD3851a	1	67165085	GRCh37	G	A	GAGGATGACA	C	T	GCAAGCTGCC		-SGIP1	ENSG00000118473	CCDS30744.1	ENST00000371037	ProteinCoding	Intronic,Promoter	r.?			None'
  def self.from_string(str)
    variant_id,
    sample_id, chr, position, genome_build,
    ref_base_plus_strand, mutant_base_plus_strand,
    five_prime_flanking_sequence_in_pyrimidine_context, ref_base_pyrimidine_context,
    mutant_base_pyrimidine_context, three_prime_flanking_sequence_in_pyrimidine_context,
    strand_of_mutation_in_pyrimidine_context,
    gene, gene_id, ccds_id, transcript_id,
    gene_type, mut_types,
    mRNA_mut_syntax, cds_mut_syntax, aa_mut_syntax,
    current_conf_status, validation_platform = str.chomp.split("\t")

    BreastCancerSNV.new(variant_id,
                        sample_id.to_sym, chr.to_sym, position.to_i, genome_build.to_sym,
                        ref_base_plus_strand.upcase.to_sym, mutant_base_plus_strand.upcase.to_sym,

                        five_prime_flanking_sequence_in_pyrimidine_context.upcase,
                        ref_base_pyrimidine_context.upcase.to_sym,
                        mutant_base_pyrimidine_context.upcase.to_sym,
                        three_prime_flanking_sequence_in_pyrimidine_context.upcase,

                        strand_of_mutation_in_pyrimidine_context.to_sym,
                        gene, gene_id, ccds_id, transcript_id,
                        gene_type.to_sym, mut_types.split(',').map(&:downcase).map(&:to_sym).to_set,
                        mRNA_mut_syntax, cds_mut_syntax, aa_mut_syntax,
                        current_conf_status.to_sym, (validation_platform || :'').to_sym )
  end


  # each_substitution_in_file('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt')
  def self.each_substitution_in_file(substitutions_filename, &block)
    return enum_for(:each_substitution_in_file, substitutions_filename).lazy  unless block_given?
    File.open(substitutions_filename) do |f|
      f.readline # skip header
      each_substitution_in_stream(f, &block)
    end
  end

  def self.each_substitution_in_stream(stream, &block)
    stream.each_line.lazy.map{|line| BreastCancerSNV.from_string(line) }.each(&block)
  end

  def mutation_types_string
    mut_types.map(&:to_s).map(&:capitalize).join(',')
  end
  private :mutation_types_string

  def to_s
    [
      variant_id, sample_id, chr, position, genome_build,
      ref_base_plus_strand, mutant_base_plus_strand,
      five_prime_flanking_sequence_in_pyrimidine_context,
      ref_base_pyrimidine_context, mutant_base_pyrimidine_context,
      three_prime_flanking_sequence_in_pyrimidine_context,
      strand_of_mutation_in_pyrimidine_context,
      gene, gene_id, ccds_id, transcript_id, gene_type,
      mutation_types_string,
      mRNA_mut_syntax, cds_mut_syntax, aa_mut_syntax,
      current_conf_status, validation_platform
    ].join("\t")
  end

  def promoter?; mut_types.include?(:promoter); end
  def intronic?; mut_types.include?(:intronic); end
  def regulatory?; promoter? || intronic?; end

  def pyrimidine_strand?
    strand_of_mutation_in_pyrimidine_context == :+
  end

  def load_sequence(genome_folder, five_prime_flank_length, three_prime_flank_length)
    File.open (File.join(genome_folder, "chr#{chr}.plain")) do |f|
      f.seek(position - five_prime_flank_length - 1)
      f.read(five_prime_flank_length + three_prime_flank_length + 1).upcase
    end
  end

  # 1-based, fully closed, given snv position is 1-based
  def interval_around_snv(five_prime_length, three_prime_flank_length)
    IntervalNotation::Syntax::Long.closed_closed(position - five_prime_length, position + three_prime_flank_length)
  end

  # 1-based, fully closed, given snv position is 1-based
  def site_interval(site, flank_length = 0)
    interval_around_snv(flank_length + site.seq_1_five_flank_length,
                        flank_length + site.seq_1_three_flank_length)
  end

  def load_site_sequence(genome_folder, site, flank_length = 0)
    load_sequence(genome_folder, flank_length + site.seq_1_five_flank_length, flank_length + site.seq_1_three_flank_length)
  end

  def five_prime_flanking_sequence_plus_strand
    pyrimidine_strand? ? five_prime_flanking_sequence_in_pyrimidine_context : revcomp(three_prime_flanking_sequence_in_pyrimidine_context)
  end

  def three_prime_flanking_sequence_plus_strand
    pyrimidine_strand? ? three_prime_flanking_sequence_in_pyrimidine_context : revcomp(five_prime_flanking_sequence_in_pyrimidine_context)
  end


  def ref_sequence_pyrimidine_context
    "#{five_prime_flanking_sequence_in_pyrimidine_context}#{ref_base_pyrimidine_context}#{three_prime_flanking_sequence_in_pyrimidine_context}"
  end

  def mutant_sequence_pyrimidine_context
    "#{five_prime_flanking_sequence_in_pyrimidine_context}#{mutant_base_pyrimidine_context}#{three_prime_flanking_sequence_in_pyrimidine_context}"
  end

  def ref_sequence_plus_strand
    "#{five_prime_flanking_sequence_plus_strand}#{ref_base_plus_strand}#{three_prime_flanking_sequence_plus_strand}"
  end

  def mutant_sequence_plus_strand
    "#{five_prime_flanking_sequence_plus_strand}#{mutant_base_plus_strand}#{three_prime_flanking_sequence_plus_strand}"
  end

  def snp_sequence_pyrimidine_context
    "#{five_prime_flanking_sequence_in_pyrimidine_context}[#{ref_base_pyrimidine_context}/#{mutant_base_pyrimidine_context}]#{three_prime_flanking_sequence_in_pyrimidine_context}"
  end

  def snp_sequence_plus_strand
    "#{five_prime_flanking_sequence_plus_strand}[#{ref_base_plus_strand}/#{mutant_base_plus_strand}]#{three_prime_flanking_sequence_plus_strand}"
  end

  def snp_sequence_from_genome(genome_folder, five_prime_flank_length, three_prime_flank_length)
    seq = load_sequence(genome_folder, five_prime_flank_length, three_prime_flank_length)

    if seq[five_prime_flank_length - five_prime_flanking_sequence_plus_strand.length, ref_sequence_plus_strand.length] == ref_sequence_plus_strand
      "#{seq[0,five_prime_flank_length]}[#{ref_base_plus_strand}/#{mutant_base_plus_strand}]#{seq[(five_prime_flank_length + 1),three_prime_flank_length]}"
    else
      raise "Error! Sequence in genome doesn't match sequence SNV flanks"
    end
  end
end

BreastCancerSNV::FILE_HEADER = [
  "variant_id", "sample_id", "chr", "position", "genome_build",
  "ref_base_plus_strand", "mutant_base_plus_strand",
  "5_prime_flanking_sequence_in_pyrimidine_context",
  "ref_base_pyrimidine_context", "mutant_base_pyrimidine_context",
  "3_prime_flanking_sequence_in_pyrimidine_context",
  "strand_of mutation_in_pyrimidine_context", "gene",
  "gene_id", "ccds_id", "transcript_id", "gene_type",
  "mut_type", "mRNA_mut_syntax", "cds_mut_syntax", "aa_mut_syntax",
  "current_conf_status", "validation_platform"
].join("\t")
