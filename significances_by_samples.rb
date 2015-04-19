$:.unshift File.absolute_path('lib', __dir__)
require 'rate_comparison_infos'
require 'statistics/statistics'

random_variants = ARGV
motif_names = File.readlines('./source_data/motif_names.txt').map(&:chomp)
sample_names = Dir.glob('./results/motif_statistics/full/any/samples/*').select{|fn| 
  File.directory?(fn) 
}.map{|fn|
  File.basename(fn)
}.reject{|fn| fn.start_with?('PD') }

sample_infos = sample_names.map{|sample_name|
  rate_infos_by_sample = random_variants.flat_map{|random_variant|
    filename = File.join('./results/motif_statistics/full/any/samples/', sample_name, "#{random_variant}.csv")
    MotifCollectionStatistics.each_in_file(filename).to_a
  }.group_by(&:motif)
  [sample_name, rate_infos_by_sample]
}.to_h


headers_top = sample_names.inject(["Motif name"]) {|sum, sample_name|
  sum + [sample_name, nil, nil, nil] 
}
puts headers_top.join("\t")

headers = sample_names.inject(["Motif name"]) {|sum, sample_name|
  sum + ["#{sample_name} cancer to random disruption ratio", "#{sample_name} cancer to random disruption ratio stddev", "#{sample_name} significance", "#{sample_name} significance stddev"] 
}
puts headers.join("\t")


motif_names.each do |motif_name|
  mean_significances = sample_names.map do |sample_name|
    infos = sample_infos[sample_name][motif_name]
    ratio_stats = Statistics.new( infos.map(&:cancer_to_random_disruption_ratio) )
    significance_stats = Statistics.new( infos.map(&:disruption_significance).map{|pvalue| -Math.log(pvalue / 0.05) } )
    [ratio_stats.mean, ratio_stats.stddev, significance_stats.mean, significance_stats.stddev]
  end

  puts ([motif_name] + mean_significances).join("\t")
end
