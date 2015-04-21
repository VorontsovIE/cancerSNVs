$:.unshift File.absolute_path('../../lib', __dir__)
require 'repeat_masker_info'
require 'data_import/breast_cancer_snv'
require 'load_genome_structure'

require 'interval_notation' # gem dependency

sites_filename = ARGV[0] # './source_data/sites_cancer.txt'
snvs_filename = ARGV[1] # './source_data/SNV_infos.txt'
genome_repeat_masker_folder = ARGV[2] # '/home/ilya/iogen/genome/hg19_repeatMasker'
raise 'Specify file with sites, SNV infos and folder with repeat masker infos'  unless sites_filename && snvs_filename && genome_repeat_masker_folder

snvs = BreastCancerSNV.each_in_file(snvs_filename).map{|snv| [snv.variant_id, snv] }.to_h
repeats_by_chromosome = read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [:Simple_repeat, :Low_complexity], expand_length: 25)
                        .map{|chromosome_name, repeats| [chromosome_name.to_s.sub(/\Achr/,'').to_sym, repeats] }.to_h

PerfectosAPE::Result.each_in_file(sites_filename).with_index.reject{|site, index|
  $stderr.puts "#{index} sites processed"  if index % 100000 == 0
  snv = snvs[site.normalized_snv_name]
  
  # it can't be correctly generalized on shuffled sequences
  # (it's hard to say whether shuffled repeats are repeats,
  # but it's more safe to remove all SNVs which have a repeat inside a 25-bp flanks)
  ### repeats_by_chromosome[snv.chromosome].intersect?(snv.site_interval(site))

  repeats_by_chromosome[snv.chromosome].include_position?(snv.position)
}.each{|site,index|
  puts site.line
}
