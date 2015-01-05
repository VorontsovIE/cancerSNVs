require 'interval_notation'
$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'
require 'cage_peak'
require 'repeat_masker_info'
require 'ensembl_exon'

def load_promoters_by_chromosome(filename, length_5_prime: 2000, length_3_prime: 500)
  promoter_intervals_by_chromosome = Hash.new{|h,k| h[k] = [] }

  File.open(filename) do |f|
    f.each_line do |line|
      chromosome, strand, tss = line.chomp.split("\t").last(3)

      chromosome_name = chromosome.to_s.start_with?('CHR') ? chromosome : "chr#{chromosome}"

      tss_pos = tss.to_i
      if strand == '1'
        promoter_intervals_by_chromosome[chromosome_name.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_5_prime, tss_pos + length_3_prime)
      elsif strand == '-1'
        promoter_intervals_by_chromosome[chromosome_name.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_3_prime, tss_pos + length_5_prime)
      end
    end
  end

  promoter_intervals_by_chromosome.map{|chromosome, intervals|
    [chromosome, IntervalNotation::Operations.union(intervals)]
  }.to_h
end

def load_cage_peaks_by_chromosome(filename, length_5_prime: 2000, length_3_prime: 500)
  CagePeak.each_in_file(filename)
          .group_by(&:chromosome)
          .map{|chromosome, peaks|
            expanded_peak_regions = peaks.map{|peak| peak.region_expanded(length_5_prime: 2000, length_3_prime: 500) }
            [chromosome, IntervalNotation::Operations.union(expanded_peak_regions)]
          }.to_h
end

def read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [], expand_length: 0)
  result = {}
  Dir.glob("#{genome_repeat_masker_folder}/*/*.fa.out") do |filename|
    chromosome_name = File.basename(filename, '.fa.out')
    repeats = RepeatMaskerInfo.each_in_file(filename).reject{|info|
      ignore_repeat_types.include?(info.repeat_type)
    }.map{|info|
      IntervalNotation::Syntax::Long.closed_closed(info.match_start - expand_length, info.match_end + expand_length)
    }.to_a # lazy iterator needs to be forced

    result[chromosome_name.to_sym] = IntervalNotation::Operations.union(repeats)
  end
  result
end

def read_coding_exons_by_chromosome(filename)
  EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              coding_regions = exons.map{|exon| exon.coding_part_region }.compact
              chromosome_name = chromosome.to_s.start_with?('chr') ? chromosome : "chr#{chromosome}"
              [chromosome_name.to_sym, IntervalNotation::Operations.union(coding_regions)]
            }.to_h
end

raise 'Specify gene TSS file' unless gene_tss_filename = ARGV[0] # './source_data/gene_tss.txt'
raise 'Specify cage peaks' unless cage_peaks_filename = ARGV[1] # '/home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt'
raise 'Specify SNV infos' unless snvs_filename = ARGV[2] # './source_data/SNV_infos_original.txt'
raise 'Specify repeat masker folder' unless genome_repeat_masker_folder = ARGV[3] # '/home/ilya/iogen/genome/hg19_repeatMasker'
raise 'Specify ensembl exons markup' unless exons_filename = ARGV[4] # '/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt'

promoters_by_chromosome = load_promoters_by_chromosome(gene_tss_filename, length_5_prime: 2000, length_3_prime: 500)
cage_peaks_by_chromosome = load_cage_peaks_by_chromosome(cage_peaks_filename, length_5_prime: 2000, length_3_prime: 500)
coding_exons_by_chromosome = read_coding_exons_by_chromosome(exons_filename)
repeats_by_chromosome = read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [:Simple_repeat, :Low_complexity], expand_length: 25)

puts BreastCancerSNV::FILE_HEADER
BreastCancerSNV.each_substitution_in_file(snvs_filename).each do |snv|
  chromosome = "chr#{snv.chr}".to_sym
  if promoters_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :promoter
  end

  if cage_peaks_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :cage_peak
  end

  if coding_exons_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :exon_coding
  end

  if repeats_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :repeat
  end

  puts snv
end
