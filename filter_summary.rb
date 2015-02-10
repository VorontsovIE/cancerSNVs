$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'optparse'
require 'rate_comparison_infos'

def load_motif_underfitting_rates(fitting_log_filename)
  File.readlines(fitting_log_filename).drop(3).slice_after("\n").map{|lines|
    motif = lines[0].chomp
    percent = lines[1].chomp.match(/\((.+)%\)$/)[1].to_f
    motif_fitting_rate = percent / 100.0
    [motif, motif_fitting_rate]
  }.to_h
end

def motifs_underfitted_in_file(fitting_log_filename, minimal_fitting_rate)
  load_motif_underfitting_rates(fitting_log_filename).select{|motif, motif_fitting_rate|
    motif_fitting_rate < minimal_fitting_rate
  }.map{|motif, motif_fitting_rate|
    motif
  }.to_set
end

motif_qualities = [:A, :B, :C]
significance_cutoff = 0.05

subjected_or_protected = :subjected
disruption_or_emergence = :disruption

minimal_fitting_rate = nil
fitting_log_filenames = []

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

  opts.on('--fitting-log FILE', 'Specify file with fitting logs in order to remove motifs',
                                'with underfitted site distributions.',
                                'One can specify several fitting logs by specifing this option several times.',
                                'Underfitting in any one log discards a motif'){|value|
    fitting_log_filenames << value
  }
  opts.on('--min-fitting-rate RATE', 'Reject motifs which were underfitted, i.e. distribution of random sites',
                                     'has less sites than a certain part of real site distribution.',
                                     'Fitting log should be turned on. Default: 0.95'){|value|
    minimal_fitting_rate = Float(value)
  }
end.parse!(ARGV)

if fitting_log_filenames.empty?
  if minimal_fitting_rate
    raise 'You specified min-fitting-rate but didn\'t specified fitting logs'
    exit 1
  end
else
  minimal_fitting_rate ||= 0.95
end

if $stdin.tty?
  raise 'Specify filename'  unless filename = ARGV[0]
  motif_infos = RateComparisonInfo.each_in_file(filename).to_a
else
  $stdin.readline  # skip header
  motif_infos = RateComparisonInfo.each_in_stream($stdin).to_a
end

if disruption_or_emergence == :disruption
  characteristic = :cancer_to_random_disruption_ratio
  characteristic_significance = :disruption_significance
else  # emergence
  characteristic = :cancer_to_random_emergence_ratio
  characteristic_significance = :emergence_significance
end

filtered_motif_infos = motif_infos.select{|infos|
  infos[characteristic]
}.select{|infos|
  if subjected_or_protected == :subjected
     infos[characteristic] > 1
  else # protected
    infos[characteristic] < 1
  end
}.select{|infos|
  infos[characteristic_significance] && infos[characteristic_significance] < significance_cutoff
}.select{|infos|
  motif_qualities.include?(infos[:quality])
}

motifs_underfitted = fitting_log_filenames.map{|fitting_log_filename|
  motifs_underfitted_in_file(fitting_log_filename, minimal_fitting_rate)
}.inject(Set.new, &:union)

filtered_motif_infos.reject!{|infos|
  motifs_underfitted.include?(infos[:motif])
}


filtered_motif_infos.sort_by!{|infos|
  infos[characteristic]
}

filtered_motif_infos.reverse!  if subjected_or_protected == :subjected

puts RateComparisonInfo::COLUMN_ORDER.map{|col| RateComparisonInfo::COLUMN_NAMES[col] }.join("\t")
filtered_motif_infos.each{|infos|
  puts infos
}
