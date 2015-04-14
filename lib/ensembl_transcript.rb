require_relative 'ensembl_exon'

class EnsemblTranscript
  attr_reader :exons
  def initialize(exons)
    raise 'No one exon specified'  if exons.empty?
    raise 'Exons from different transcripts'  unless exons.map(&:ensembl_transcript_id).uniq.size == 1
    raise 'Exons from different chromosomes'  unless exons.map(&:chromosome).uniq.size == 1
    raise 'Exons from different strands'  unless exons.map(&:strand).uniq.size == 1
    raise 'Exons have different gene_biotype'  unless exons.map(&:gene_biotype).uniq.size == 1
    raise 'Exons have different associated_gene_name'  unless exons.map(&:associated_gene_name).uniq.size == 1
    raise 'Exons have different cds_length'  unless exons.map(&:cds_length).uniq.size == 1
    @exons = exons
  end
  def ensembl_transcript_id
    @exons.first.ensembl_transcript_id
  end
  def ensembl_gene_id
    @exons.first.ensembl_gene_id
  end
  def strand
    @exons.first.strand
  end
  def chromosome
    @exons.first.chromosome
  end
  def gene_biotype
    @exons.first.gene_biotype
  end
  def associated_gene_name
    @exons.first.associated_gene_name
  end
  def description
    @exons.first.description
  end
  def exonic_region
    IntervalNotation::Operations.union( @exons.map(&:exon_region) )
  end

  def cds_length
    @exons.first.cds_length
  end

  def has_coding_part?
    !!(cds_length && cds_length > 0)
  end

  def has_coding_exons?
    @exons.any?{|exon| exon.phase != -1 }
  end

  # Only exonic part of UTR
  def utr_5
    IntervalNotation::Operations.union( @exons.map(&:utr_5) )
  end
  def utr_3
    IntervalNotation::Operations.union( @exons.map(&:utr_3) )
  end
  def coding_region
    IntervalNotation::Operations.union( @exons.map(&:coding_part_region) )
  end
  def exonic_non_utr
    exonic_region - utr_3 - utr_5
  end
  def transcript_start
    exons.first.transcript_start # any exon has the same transcript start
  end
  def transcript_end
    exons.first.transcript_end
  end
  def promoter(five_prime_length:, three_prime_length:)
    case strand
    when :+
      IntervalNotation::Syntax::Long.closed_closed(transcript_start - five_prime_length, transcript_start + three_prime_length - 1)
    when :-
      IntervalNotation::Syntax::Long.closed_closed(transcript_end - three_prime_length + 1, transcript_end + five_prime_length)
    else
      raise 'Wrong strand format'
    end
  end
end
