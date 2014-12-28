$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'
require 'interval_notation'

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

gene_tss_filename = ARGV[0] # 'source_data/gene_tss.txt'
snvs_filename = ARGV[1] #'source_data/SUBSTITUTIONS_13Apr2012_snz.txt'

raise 'Specify file with gene TSSes and file with SNV infos'  unless gene_tss_filename && snvs_filename

promoters_by_chromosome = load_promoters_by_chromosome(gene_tss_filename, length_5_prime: 2000, length_3_prime: 500)

puts BreastCancerSNV::FILE_HEADER
BreastCancerSNV.each_substitution_in_file(snvs_filename).each do |snv|
  if promoters_by_chromosome[snv.chr.to_sym].include_position?(snv.position)
    snv.mut_types << :promoter
  end
  puts snv
end
