def to_float(str)
  str && str.to_f
end

lines = ARGF.readlines
puts lines.shift

lines.map{|line|
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
}.select{|infos|
  infos[:disruption_significance] && infos[:disruption_significance] < 0.05 || 
  infos[:emergence_significance] && infos[:emergence_significance] < 0.05
}.reject{|infos| 
 # infos[:quality] == :D
  false
}.sort_by{|infos| infos[:cancer_to_random_disruption_ratio] }.reverse.each{|infos|
  puts infos[:line]
}
