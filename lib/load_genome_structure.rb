require 'interval_notation'
require_relative 'cage_peak'
require_relative 'ensembl_exon'

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def load_promoters_by_chromosome(filename, length_5_prime: 5000, length_3_prime: 500)
  result = EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              promoter_regions = exons.select(&:transcript_start).map{|exon|
                tss_pos = exon.transcript_start
                if exon.strand == :+
                  IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_5_prime, tss_pos + length_3_prime)
                else
                  IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_3_prime, tss_pos + length_5_prime)
                end
              }

              [chromosome, IntervalNotation::Operations.union(promoter_regions)]
            }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end

def load_cage_peaks_by_chromosome(filename, length_5_prime: 2000, length_3_prime: 500)
  CagePeak.each_in_file(filename)
          .group_by(&:chromosome)
          .map{|chromosome, peaks|
            expanded_peak_regions = peaks.map{|peak| peak.region_expanded(length_5_prime: length_5_prime, length_3_prime: length_3_prime) }
            [chromosome, IntervalNotation::Operations.union(expanded_peak_regions)]
          }.to_h
end

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def read_coding_exons_by_chromosome(filename)
  EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              [chromosome, IntervalNotation::Operations.union(exons.map(&:coding_part_region))]
            }.to_h
end

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def read_introns_by_chromosome(filename)
  result = EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              transcript_introns = exons.group_by(&:ensembl_transcript_id).map{|ensembl_transcript_id, exons|
                gene_exons = IntervalNotation::Operations.union(exons.map(&:exon_region))
                gene_exons.covering_interval - gene_exons
              }
              [chromosome, IntervalNotation::Operations.union(transcript_introns)]
            }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end

### Regions of kataegis
KataegisRegion = Struct.new(:cancer_type, :sample_name, :chromosome, :position_start, :position_end, :number_of_variants) do
  def interval(expansion_length: 0)
    @interval ||= IntervalNotation::Syntax::Long.closed_closed(position_start - expansion_length, position_end + expansion_length)
  end

  def to_s
    [cancer_type, sample_name, chromosome, position_start, position_end, number_of_variants].join("\t")
  end

  # Breast Cancer PD4224a 1 8976832 8982950 9
  def self.from_string(line)
    cancer_type, sample_name, chromosome, position_start, position_end, number_of_variants = line.chomp.split("\t")
    self.new( cancer_type.to_sym, sample_name.to_sym,
              chromosome.to_sym, position_start.to_i, position_end.to_i,
              number_of_variants.to_i )
  end

  def self.each_in_file(filename, &block)
    File.open(filename){|f|
      f.readline
      f.each_line.map{|line| self.from_string(line) }.each(&block)
    }
  end
end

# ./source_data/AlexandrovEtAl/coordinates_of_kataegis.csv
def load_kataegis_regions_by_chromosome(filename, expansion_length: 1000)
  kataegis_regions = KataegisRegion.each_in_file(filename).to_a
  result = kataegis_regions.group_by(&:chromosome).map{|chromosome, regions|
    expanded_intervals = regions.map{|region| region.interval(expansion_length: expansion_length) }

    [chromosome, IntervalNotation::Operations.union(expanded_intervals)]
  }.to_h

  result.default = IntervalNotation::Syntax::Long::Empty
  result
end
