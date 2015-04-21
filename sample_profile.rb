$:.unshift File.absolute_path('lib', __dir__)
require 'rate_comparison_infos'
require 'statistics/statistics'
require '../experiment_configuration'
require 'shellwords'
require 'optparse'

name = nil
OptionParser.new{|opts|
  opts.banner = "Usage: #{opts.program_name} <list of motif statistics files> [options]\n" +
                "         or\n" +
                "       <list of motif statistics files> | #{opts.program_name} [options]"
  opts.on('--sample-name NAME', 'Specify sample name. It will be used as a header'){|value| name = value }
}.parse!(ARGV)

class Array
  def mean
    self.inject(0.0, &:+) / self.size
  end
  def variance
    m = self.mean
    self.map{|x| (x-m)**2 }.inject(0.0, &:+) / self.size
  end
  def stddev
    Math.sqrt(self.variance * self.size / (self.size - 1))
  end
end

if $stdin.tty?
  filelist = ARGV
else
  filelist = $stdin.read.shellsplit
end

all_genome_motif_infos = filelist.flat_map do |filename|
  MotifCollectionStatistics.each_in_file(filename).to_a#.map{|info| [info.motif, info] }.to_h
end.group_by(&:motif)

motif_names = File.readlines(MOTIF_NAMES_FILE).map(&:chomp)

# motif_names.each do |motif|
#   disruption_rates = all_genome_motif_infos[motif].select{|infos|
#     infos.disruption_significance <= 0.05
#   }.map(&:cancer_to_random_disruption_ratio)
  
#   stat = Statistics.new(disruption_rates)
#   if stat.relative_stddev #&& stat.relative_stddev > 0.2
#     puts [motif,
#           stat.values.map{|x| x.round(3)}.inspect,
#           stat.mean && stat.mean.round(3),
#           stat.stddev && stat.stddev.round(3),
#           stat.relative_stddev ? "#{(100 * stat.relative_stddev).round(2)}%" : nil].join("\t")
#   end
# end

sample_average_profile = motif_names.map do |motif|
  motif_disruption_rates = all_genome_motif_infos[motif].map(&:cancer_to_random_disruption_ratio)
  rate_stat = Statistics.new(motif_disruption_rates)
  rate_stat.mean
end

sample_average_weights = motif_names.map do |motif|
  motif_disruption_rates = all_genome_motif_infos[motif].map(&:cancer_to_random_disruption_ratio)
  motif_weights = all_genome_motif_infos[motif].map(&:disruption_significance).map{|el| - Math.log(el) }
  weights_stat = Statistics.new(motif_weights)
  weights_stat.mean
end

puts [">#{name}", *sample_average_weights].join("\t")
puts [name, *sample_average_profile].join("\t")
