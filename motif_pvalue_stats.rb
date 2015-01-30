require 'shellwords'

def to_float(str)
  str && str.to_f
end

def mean(values)
  values.inject(0.0, &:+) / values.size
end

def geom_mean(values)
  values.inject(1.0, &:*) ** (1.0/ values.size)
end

def variance(values)
  m = mean(values)
  values.map{|x| (x - m) ** 2 }.inject(0.0, &:+)
end

def stddev(values)
  Math.sqrt(variance(values) / (values.size - 1))
end

if $stdin.tty?
  files = ARGV
else
  files = Shellwords.split($stdin.read)
end

motif_infos = files.flat_map do |fn|
  File.readlines(fn).drop(1).map do |line|
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
  end
end

puts ['Motif',
      'Disruption P-value mean', 'Disruption P-value mean (geometric)', 'Disruption P-value stddev',
      'Emergence P-value mean', 'Emergence P-value mean (geometric)', 'Emergence P-value stddev'
    ].join("\t")
motif_infos.group_by{|info| info[:motif] }.sort.each do |motif, infos|
  disruption_pvalues = infos.map{|info| info[:disruption_significance] }
  emergence_pvalues = infos.map{|info| info[:emergence_significance] }
  puts [motif,
        mean(disruption_pvalues), geom_mean(disruption_pvalues), stddev(disruption_pvalues),
        mean(emergence_pvalues), geom_mean(emergence_pvalues), stddev(emergence_pvalues)
      ].join("\t")
end
