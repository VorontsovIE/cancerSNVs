$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'optparse'
require 'rate_comparison_infos'

motif_qualities = [:A, :B, :C]
significance_cutoff = 0.05

subjected_or_protected = :subjected
disruption_or_emergence = :disruption

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <file to filter> [options]"
  opts.separator('Options:')
  opts.on('--motif-qualities QUALITIES', 'Take motif of given qualities (comma-separated). Default is A,B,C') {|value|
    motif_qualities = value.upcase.split(',').map(&:to_sym)
  }

  opts.on('--subjected', 'Select only sites subjected to site disruption/emergence.', 'Don\'t use with --protected'){ subjected_or_protected = :subjected }
  opts.on('--protected', 'Select only sites protected from site disruption/emergence.', 'Don\'t use with --subjected'){ subjected_or_protected = :protected }

  opts.on('--disruption', 'Select only sites subjected or protected to site disruption.', 'Don\'t use with --emergence'){ disruption_or_emergence = :disruption }
  opts.on('--emergence', 'Select only sites subjected or protected to site emergence.', 'Don\'t use with --disruption'){ disruption_or_emergence = :emergence }

  opts.on('--significance CUTOFF', 'Significance cutoff. Default is 0.05') {|value|
    significance_cutoff = Float(value)
  }
end.parse!(ARGV)

if $stdin.tty?
  raise 'Specify filename'  unless filename = ARGV[0]
  motif_infos = MotifCollectionStatistics.each_in_file(filename).to_a
else
  $stdin.readline  # skip header
  motif_infos = MotifCollectionStatistics.each_in_stream($stdin).to_a
end

if disruption_or_emergence == :disruption
  characteristic = :cancer_to_random_disruption_ratio
  characteristic_significance = :disruption_significance_fitting_aware
else  # emergence
  characteristic = :cancer_to_random_emergence_ratio
  characteristic_significance = :emergence_significance_fitting_aware
end

filtered_motif_infos = motif_infos.select{|infos|
  infos.send(characteristic)
}.select{|infos|
  if subjected_or_protected == :subjected
     infos.send(characteristic) > 1
  else # protected
    infos.send(characteristic) < 1
  end
}.select{|infos|
  infos.send(characteristic_significance) && infos.send(characteristic_significance) < significance_cutoff
}.select{|infos|
  motif_qualities.include?(infos.quality)
}

filtered_motif_infos.sort_by!{|infos|
  infos.send(characteristic)
}

filtered_motif_infos.reverse!  if subjected_or_protected == :subjected

puts MotifCollectionStatistics.table_header
filtered_motif_infos.each{|infos|
  puts infos
}
