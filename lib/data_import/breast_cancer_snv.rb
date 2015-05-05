require 'forwardable'
require 'interval_notation'
require_relative '../region_type'
require_relative '../sequence'
require_relative '../sequence_with_snv'

BreastCancerSNV = Struct.new( :variant_id,
                              :sample_id, :chromosome, :position, :genome_build,

                              :ref_base_plus_strand, :mutant_base_plus_strand,

                              :five_prime_flanking_sequence_in_pyrimidine_context,
                              :ref_base_pyrimidine_context,
                              :mutant_base_pyrimidine_context,
                              :three_prime_flanking_sequence_in_pyrimidine_context,

                              :strand_of_mutation_in_pyrimidine_context,
                              :gene, :gene_id, :ccds_id, :transcript_id,

                              :gene_type, :mutation_region_types,

                              :mRNA_mut_syntax, :cds_mut_syntax, :aa_mut_syntax,
                              :current_conf_status, :validation_platform ) do


  attr_accessor :_mutation_region_types_raw

  extend Forwardable
  def_delegators :mutation_region_types, *RegionType::FEATURE_INQUIRIES

  def mutation_region_types
    @mutation_region_types ||= RegionType.from_string(_mutation_region_types_raw)
  end

  # '27607140	PD3851a	1	67165085	GRCh37	G	A	GAGGATGACA	C	T	GCAAGCTGCC		-SGIP1	ENSG00000118473	CCDS30744.1	ENST00000371037	ProteinCoding	Intronic,Promoter	r.?			None'
  def self.from_string(str)
    variant_id,
    sample_id, chromosome, position, genome_build,
    ref_base_plus_strand, mutant_base_plus_strand,
    five_prime_flanking_sequence_in_pyrimidine_context, ref_base_pyrimidine_context,
    mutant_base_pyrimidine_context, three_prime_flanking_sequence_in_pyrimidine_context,
    strand_of_mutation_in_pyrimidine_context,
    gene, gene_id, ccds_id, transcript_id,
    gene_type, mutation_region_types_raw,
    mRNA_mut_syntax, cds_mut_syntax, aa_mut_syntax,
    current_conf_status, validation_platform = str.chomp.split("\t", 23)

    BreastCancerSNV.new(variant_id,
                        sample_id.to_sym, chromosome.to_sym, position.to_i, genome_build.to_sym,
                        ref_base_plus_strand.upcase.to_sym, mutant_base_plus_strand.upcase.to_sym,

                        five_prime_flanking_sequence_in_pyrimidine_context.upcase,
                        ref_base_pyrimidine_context.upcase.to_sym,
                        mutant_base_pyrimidine_context.upcase.to_sym,
                        three_prime_flanking_sequence_in_pyrimidine_context.upcase,

                        strand_of_mutation_in_pyrimidine_context.to_sym,
                        gene, gene_id, ccds_id, transcript_id,
                        gene_type.to_sym, nil, # don't initialize mut_type right now
                        mRNA_mut_syntax, cds_mut_syntax, aa_mut_syntax,
                        current_conf_status.to_sym, validation_platform.to_sym )
    .tap {|snv|
      snv._mutation_region_types_raw = mutation_region_types_raw # lazy initialization of mutation_region_types
    }
  end


  # each_in_file('./source_data/SNV_infos.txt')
  def self.each_in_file(substitutions_filename, &block)
    return enum_for(:each_in_file, substitutions_filename).lazy  unless block_given?
    File.open(substitutions_filename) do |f|
      f.readline # skip header
      each_in_stream(f, &block)
    end
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.lazy.map{|line| self.from_string(line) }.each(&block)
  end

  def to_s
    [
      variant_id, sample_id, chromosome, position, genome_build,
      ref_base_plus_strand, mutant_base_plus_strand,
      five_prime_flanking_sequence_in_pyrimidine_context,
      ref_base_pyrimidine_context, mutant_base_pyrimidine_context,
      three_prime_flanking_sequence_in_pyrimidine_context,
      strand_of_mutation_in_pyrimidine_context,
      gene, gene_id, ccds_id, transcript_id, gene_type,
      mutation_region_types,
      mRNA_mut_syntax, cds_mut_syntax, aa_mut_syntax,
      current_conf_status, validation_platform
    ].join("\t")
  end

  def pyrimidine_strand?
    strand_of_mutation_in_pyrimidine_context == :+
  end

  def load_sequence(genome_reader, five_prime_flank_length, three_prime_flank_length)
    genome_reader.load_sequence(chromosome, ONE_BASED_INCLUSIVE, position - five_prime_flank_length, position + three_prime_flank_length)
  end

  # 1-based, fully closed, given snv position is 1-based
  def interval_around_snv(five_prime_length, three_prime_flank_length)
    IntervalNotation::Syntax::Long.closed_closed(position - five_prime_length, position + three_prime_flank_length)
  end

  # Deprecated
  # 1-based, fully closed, given snv position is 1-based
  def site_interval(site, flank_length = 0)
    interval_around_snv(flank_length + site.seq_1_five_flank_length,
                        flank_length + site.seq_1_three_flank_length)
  end

  # Deprecated
  def load_site_sequence(genome_reader, site, flank_length = 0)
    load_sequence(genome_reader, flank_length + site.seq_1_five_flank_length, flank_length + site.seq_1_three_flank_length)
  end

  def five_prime_flanking_sequence_plus_strand
    pyrimidine_strand? ? five_prime_flanking_sequence_in_pyrimidine_context : Sequence.revcomp(three_prime_flanking_sequence_in_pyrimidine_context)
  end

  def three_prime_flanking_sequence_plus_strand
    pyrimidine_strand? ? three_prime_flanking_sequence_in_pyrimidine_context : Sequence.revcomp(five_prime_flanking_sequence_in_pyrimidine_context)
  end

  def snv_sequence_pyrimidine_context(five_prime_flank_length: nil, three_prime_flank_length: nil)
    five_prime_flank_length ||= five_prime_flanking_sequence_in_pyrimidine_context.length
    three_prime_flank_length ||= three_prime_flanking_sequence_in_pyrimidine_context.length
    five_prime_chunk = five_prime_flanking_sequence_in_pyrimidine_context[-five_prime_flank_length..-1]
    three_prime_chunk = three_prime_flanking_sequence_in_pyrimidine_context[0, three_prime_flank_length]
    "#{ five_prime_chunk }[#{ref_base_pyrimidine_context}/#{mutant_base_pyrimidine_context}]#{ three_prime_chunk }"
  end

  # SNV sequence from information in object itself (doesn't involve loading sequence from genome). See also #snv_sequence_from_genome
  def sequence_with_snv(five_prime_flank_length: nil, three_prime_flank_length: nil)
    five_prime_flank_length ||= five_prime_flanking_sequence_plus_strand.length
    three_prime_flank_length ||= three_prime_flanking_sequence_plus_strand.length
    if five_prime_flank_length > five_prime_flanking_sequence_plus_strand.length
      raise "Can't obtain such a long subsequence from provided data: requested #{five_prime_flank_length}bp. But provided 5' flank is of length #{five_prime_flanking_sequence_plus_strand.length}.\n" +
            "Use snv_sequence_from_genome instead"
    end
    if three_prime_flank_length > three_prime_flanking_sequence_plus_strand.length
      raise "Can't obtain such a long subsequence from provided data: requested #{three_prime_flank_length}bp. But provided 3' flank is of length #{three_prime_flanking_sequence_plus_strand.length}.\n" +
            "Use snv_sequence_from_genome instead"
    end
    seq_5 = five_prime_flanking_sequence_plus_strand[five_prime_flanking_sequence_plus_strand.length - five_prime_flank_length, five_prime_flank_length]
    seq_3 = three_prime_flanking_sequence_plus_strand[0, three_prime_flank_length]
    allele_variants = [ref_base_plus_strand, mutant_base_plus_strand]
    SequenceWithSNV.new(seq_5, allele_variants, seq_3)
  end

  # Not used at a moment
  def snv_sequence_from_genome(genome_reader, five_prime_flank_length, three_prime_flank_length)
    seq = load_sequence(genome_reader, five_prime_flank_length, three_prime_flank_length)

    seq_5 = seq[0, five_prime_flank_length]
    seq_3 = seq[(five_prime_flank_length + 1), three_prime_flank_length]
    allele_variants = [ref_base_plus_strand, mutant_base_plus_strand]
    SequenceWithSNV.new(seq_5, allele_variants, seq_3)
  end

  # 1-bp context on plus strand
  def context_before_snv_plus_strand
    @context_before_snv_plus_strand ||= "#{five_prime_flanking_sequence_plus_strand[-1]}#{ref_base_plus_strand}#{three_prime_flanking_sequence_plus_strand[0]}"
  end

  def to_snv_info(genome_reader, flank_length: 25)
    seq = genome_reader.read_sequence(chromosome, ONE_BASED_INCLUSIVE, position - flank_length, position + flank_length).upcase
    before_substitution = ref_base_plus_strand
    after_substitution = mutant_base_plus_strand
    unless seq[flank_length] == before_substitution.to_s
      raise "Base before substitution `#{before_substitution}` is not consistent with reference genome `#{seq[flank_length]}`"
    end
    five_prime_flank = seq[0, flank_length]
    three_prime_flank = seq[flank_length + 1, flank_length]
    allele_variants = [before_substitution, after_substitution]
    snv_seq = SequenceWithSNV.new(five_prime_flank, allele_variants, three_prime_flank)

    SNVInfo.new(variant_id, snv_seq, :'Breast cancer (Nik-Zainal)', sample_id, chromosome, position, :+, mutation_region_types)
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
