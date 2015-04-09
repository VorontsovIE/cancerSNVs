$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'optparse'
require 'rate_comparison_infos'

motif_qualities = [:A, :B, :C]
significance_cutoff = 0.05

subjected_or_protected = :subjected
characteristic = nil
characteristic_significance = nil

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <file to filter> [options]"
  opts.separator('Options:')
  opts.on('--motif-qualities QUALITIES', 'Take motifs of given qualities (comma-separated). Default is A,B,C') {|value|
    motif_qualities = value.upcase.split(',').map(&:to_sym)
  }

  opts.on('--subjected', 'Select only sites subjected to events specified in characteristic. This is default.', 'Don\'t use with --protected'){ subjected_or_protected = :subjected }
  opts.on('--protected', 'Select only sites protected to events specified in characteristic', 'Don\'t use with --subjected'){ subjected_or_protected = :protected }

  opts.on('--characteristic MODE', 'MODE can be disruption/emergence/substitution-in-core. Select only sites subjected or protected to this kind of events') do |mode|
    case mode.downcase
    when 'disruption'
      characteristic = :cancer_to_random_disruption_ratio
      characteristic_significance = :disruption_significance_fitting_aware
    when 'emergence'
      characteristic = :cancer_to_random_emergence_ratio
      characteristic_significance = :emergence_significance_fitting_aware
    when 'substitution-in-core'
      characteristic = :cancer_to_random_core_ratio
      characteristic_significance = :core_flank_significance_fitting_aware
    else
      raise 'Unknown characteristic'
    end
  end

  opts.on('--significance CUTOFF', 'Significance cutoff for characteristic significance. Default is 0.05') {|value|
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

if characteristic
  filtered_motif_infos = motif_infos.select{|infos|
    infos.send(characteristic)
  }.select{|infos|
    case subjected_or_protected
    when :subjected
      infos.send(characteristic) > 1
    when :protected
      infos.send(characteristic) < 1
    else
      true
    end
  }.select{|infos|
    infos.send(characteristic_significance) && infos.send(characteristic_significance) < significance_cutoff
  }
else
  filtered_motif_infos = motif_infos
end

filtered_motif_infos = filtered_motif_infos.select{|infos|
  motif_qualities.include?(infos.quality)
}

if characteristic
  filtered_motif_infos.sort_by!{|infos|
    infos.send(characteristic)
  }
  filtered_motif_infos.reverse!  if subjected_or_protected == :subjected
else
  filtered_motif_infos.sort_by!{|infos| infos.motif }
end


puts MotifCollectionStatistics.table_header
filtered_motif_infos.each{|infos|
  puts infos
}
