$:.unshift File.absolute_path('../../lib', __dir__)
require 'data_import/breast_cancer_snv'
require 'load_genome_structure'

raise 'Specify SNV infos' unless snvs_filename = ARGV[0] # './source_data/SNV_infos_original.txt'
raise 'Specify ensembl exons markup' unless exons_filename = ARGV[1] # '/home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt'

promoters_by_chromosome = load_promoters_by_chromosome(exons_filename, length_5_prime: 2000, length_3_prime: 500)
introns_by_chromosome = read_introns_by_chromosome(exons_filename)

puts BreastCancerSNV::FILE_HEADER
BreastCancerSNV.each_in_file(snvs_filename).each do |snv|
  chromosome = "chr#{snv.chromosome}".to_sym
  snv.mutation_region_types.discard_types! # discard all mutation types
  if promoters_by_chromosome[chromosome].include_position?(snv.position)
    snv.mutation_region_types << :promoter
  end

  if introns_by_chromosome[chromosome].include_position?(snv.position)
    snv.mutation_region_types << :intronic
  end

  puts snv
end
