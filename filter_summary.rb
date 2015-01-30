require 'optparse'

def to_float(str)
  str && str.to_f
end

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
    significance_cutoff = value.to_f
  }
end.parse!(ARGV)

lines = ARGF.readlines

puts lines.shift

motif_infos = lines.map{|line|
  motif, official_gene_name, quality,
          cancer_to_random_disruption_ratio, disruption_significance,
          cancer_to_random_emergence_ratio, emergence_significance,
          random_disrupted, random_emerged, random_total_before_substitution, random_total_after_substitution,
          cancer_disrupted, cancer_emerged, cancer_total_before_substitution, cancer_total_after_substitution,
          cancer_disruption_rate, cancer_emergence_rate, random_disruption_rate, random_emergence_rate,
          disruption_significance_uncorrected, emergence_significance_uncorrected,
          gene = line.chomp.split("\t")

  {
    line: line, motif: motif,
    official_gene_name: official_gene_name, quality: quality.upcase.to_sym,
    cancer_to_random_disruption_ratio: to_float(cancer_to_random_disruption_ratio), disruption_significance: to_float(disruption_significance),
    cancer_to_random_emergence_ratio: to_float(cancer_to_random_emergence_ratio), emergence_significance: to_float(emergence_significance),
    random_disrupted: random_disrupted.to_i, random_emerged: random_emerged.to_i,
    random_total_before_substitution: random_total_before_substitution.to_i, random_total_after_substitution: random_total_after_substitution.to_i,
    cancer_disrupted: cancer_disrupted.to_i, cancer_emerged: cancer_emerged.to_i,
    cancer_total_before_substitution: cancer_total_before_substitution.to_i, cancer_total_after_substitution: cancer_total_after_substitution.to_i,
    cancer_disruption_rate: cancer_disruption_rate.to_f, cancer_emergence_rate: cancer_emergence_rate.to_f,
    random_disruption_rate: random_disruption_rate.to_f, random_emergence_rate: random_emergence_rate.to_f,
    disruption_significance_uncorrected: to_float(disruption_significance_uncorrected),
    emergence_significance_uncorrected: to_float(emergence_significance_uncorrected),
    gene: gene,
  }
}

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
}.sort_by{|infos|
  infos[characteristic]
}

filtered_motif_infos.reverse!  if subjected_or_protected == :subjected

filtered_motif_infos.each{|infos|
  puts infos[:line]
}
