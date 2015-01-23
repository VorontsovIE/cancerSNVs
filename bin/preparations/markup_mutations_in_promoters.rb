require 'interval_notation'
$:.unshift File.absolute_path('../../lib', __dir__)
require 'breast_cancer_snv'
require 'cage_peak'
require 'repeat_masker_info'
require 'ensembl_exon'
require 'load_genome_structure'

raise 'Specify gene TSS file' unless gene_tss_filename = ARGV[0] # './source_data/gene_tss.txt'
raise 'Specify cage peaks' unless cage_peaks_filename = ARGV[1] # '/home/ilya/iogen/cages/hg19/freeze1/hg19.cage_peak_tpm_ann.osc.txt'
raise 'Specify SNV infos' unless snvs_filename = ARGV[2] # './source_data/SNV_infos_original.txt'
raise 'Specify repeat masker folder' unless genome_repeat_masker_folder = ARGV[3] # '/home/ilya/iogen/genome/hg19_repeatMasker'
raise 'Specify ensembl exons markup' unless exons_filename = ARGV[4] # '/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt'

promoters_by_chromosome = load_promoters_by_chromosome(gene_tss_filename, length_5_prime: 2000, length_3_prime: 500)
cage_peaks_by_chromosome = load_cage_peaks_by_chromosome(cage_peaks_filename, length_5_prime: 2000, length_3_prime: 500)
coding_exons_by_chromosome = read_coding_exons_by_chromosome(exons_filename)
introns_by_chromosome = read_introns_by_chromosome(exons_filename)
repeats_by_chromosome = read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [:Simple_repeat, :Low_complexity], expand_length: 25)

puts BreastCancerSNV::FILE_HEADER
BreastCancerSNV.each_substitution_in_file(snvs_filename).each do |snv|
  chromosome = "chr#{snv.chr}".to_sym
  if promoters_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :promoter
  end

  if introns_by_chromosome[chromosome].include_position?(snv.position)
    snv.mut_types << :intronic
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
