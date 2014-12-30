$:.unshift File.absolute_path('../../lib', __dir__)
require 'repeat_masker_info'
require 'breast_cancer_snv'

require 'interval_notation' # gem dependency

def read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [])
  result = {}
  Dir.glob("#{genome_repeat_masker_folder}/*/*.fa.out") do |filename|
    chromosome_name = File.basename(filename, '.fa.out')
    $stderr.puts chromosome_name
    repeats = RepeatMaskerInfo.each_in_file(filename).reject{|info|
      ignore_repeat_types.include?(info.repeat_type)
    }.map{|info|
      IntervalNotation::Syntax::Long.closed_closed(info.match_start, info.match_end)
    }.to_a # lazy iterator needs to be forced
    result[chromosome_name.sub(/\Achr/,'').to_sym] = IntervalNotation::Operations.union(repeats)  
  end
  result
end

snvs_filename = 'source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt'
genome_repeat_masker_folder = '/home/ilya/iogen/genome/hg19_repeatMasker'
sites_filename = ARGV.first # 'source_data/sites_cancer.txt'
raise 'Specify file with sites'  unless sites_filename

tm = Time.now
snvs = BreastCancerSNV.each_substitution_in_file(snvs_filename).map{|snv| [snv.variant_id, snv] }.to_h
$stderr.puts("SNV infos loaded in #{Time.now - tm} sec.")
tm = Time.now
repeats_by_chromosome = read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [:Simple_repeat, :Low_complexity])
$stderr.puts("repeats loaded in #{Time.now - tm} sec.")


MutatatedSiteInfo.each_site(sites_filename).with_index.reject{|site, index|
  $stderr.puts "#{index} sites processed"  if index % 100000 == 0
  snv = snvs[site.normalized_snp_name]
  
  # it can't be correctly generalized on shuffled sequences
  # (it's hard to say whether shuffled repeats are repeats,
  # but it's more safe to remove all SNVs which have a repeat inside a 25-bp flanks)
  ### repeats_by_chromosome[snv.chr].intersect?(snv.site_interval(site))

  repeats_by_chromosome[snv.chr].intersect?(snv.interval_around_snv(25,25))
}.each{|site,index|
  puts site.line
}
