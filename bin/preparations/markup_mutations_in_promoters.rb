$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'
require 'interval_notation'
require 'cage_peak'

def load_promoters_by_chromosome(filename, length_5_prime: 2000, length_3_prime: 500)
  promoter_intervals_by_chromosome = Hash.new{|h,k| h[k] = [] }

  File.open(filename) do |f|
    f.each_line do |line|
      chr, strand, tss = line.chomp.split("\t").last(3)
      tss_pos = tss.to_i
      if strand == '1'
        promoter_intervals_by_chromosome[chr.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_5_prime, tss_pos + length_3_prime)
      elsif strand == '-1'
        promoter_intervals_by_chromosome[chr.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_3_prime, tss_pos + length_5_prime)
      end
    end
  end

  promoter_intervals_by_chromosome.map{|chr, intervals|
    [chr, IntervalNotation::Operations.union(intervals)]
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

gene_tss_filename = ARGV[0] # './source_data/gene_tss.txt'
cage_peaks_filename = ARGV[1] # '/home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt'
snvs_filename = ARGV[2] # './source_data/SNV_infos_original.txt'

raise 'Specify file with gene TSSes and file with SNV infos'  unless gene_tss_filename && snvs_filename

promoters_by_chromosome = load_promoters_by_chromosome(gene_tss_filename, length_5_prime: 2000, length_3_prime: 500)
cage_peaks_by_chromosome = load_cage_peaks_by_chromosome(cage_peaks_filename, length_5_prime: 2000, length_3_prime: 500)
# TODO: We should cut coding exons part from expanded cage peaks!!!

puts BreastCancerSNV::FILE_HEADER
BreastCancerSNV.each_substitution_in_file(snvs_filename).each do |snv|
  chromosome = snv.chr.to_sym
  if promoters_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :promoter
  end

  if cage_peaks_by_chromosome["chr#{chromosome}".to_sym].include_position?(snv.position)
    snv.mut_types << :cage_peak
  end
  puts snv
end
