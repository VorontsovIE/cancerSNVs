def to_float(str)
  str && str.to_f
end

# # unclassified is not included in total
# FischerTableData = Struct.new(:first_total, :first_positive, :second_total, :second_positive, :first_unclassified, :second_unclassified) do
#   def first_negative
#     first_total - first_positive
#   end
#   def second_negative
#     second_total - second_positive
#   end

#   def first_rate
#     first_positive.to_f / first_total
#   end
#   def second_rate
#     second_positive.to_f / second_total
#   end
#   def first_to_second_ratio
#     first_rate / second_rate
#   end

#   def significance # aware of unclassified items
#     PvalueCalculator
#       .new(class_counts: :class_and_total)
#       .calculate(
#         [first_positive, second_positive, first_total, second_total],
#         unclassified_a: first_unclassified,
#         unclassified_b: second_unclassified
#       )
#   end
# end

# RateComparisonInfo_2 = Struct.new(
#   :motif,
#   :disruption_rates,
#   :emergence_rates,
#   :official_gene_name, :quality, :gene
# ) do

#   def self.from_string(line)
#     motif, official_gene_name, quality,
#     cancer_to_random_disruption_ratio, disruption_significance_fitting_aware,
#     cancer_to_random_emergence_ratio, emergence_significance_fitting_aware,
#     random_disrupted, random_emerged, random_total_before_substitution, random_total_after_substitution,
#     cancer_disrupted, cancer_emerged, cancer_total_before_substitution, cancer_total_after_substitution,
#     cancer_disruption_rate, cancer_emergence_rate, random_disruption_rate, random_emergence_rate,
#     disruption_significance_uncorrected, emergence_significance_uncorrected,
#     gene, underfitted,
#     disruption_significance, emergence_significance,
#     disruption_significance_uncorrected_fitting_aware, emergence_significance_uncorrected_fitting_aware = line.chomp.split("\t")
#     RateComparisonInfo_2.new(
#         motif,
#         FischerTableData.new(cancer_total_before.to_i, cancer_disrupted.to_i, random_total_before.to_i, random_disrupted.to_i, 0, underfitted.to_i)
#         FischerTableData.new(cancer_total_after.to_i,  cancer_emerged.to_i,   random_total_after.to_i,  random_emerged.to_i,   0, underfitted.to_i)
#         official_gene_name, quality.upcase.to_sym, gene
#     )
#   end
# end

RateComparisonInfo = Struct.new(
                  :motif, :official_gene_name, :quality,
                  :cancer_to_random_disruption_ratio, :disruption_significance_fitting_aware,
                  :cancer_to_random_emergence_ratio, :emergence_significance_fitting_aware,
                  :random_disrupted, :random_emerged, :random_total_before_substitution, :random_total_after_substitution,
                  :cancer_disrupted, :cancer_emerged, :cancer_total_before_substitution, :cancer_total_after_substitution,
                  :cancer_disruption_rate, :cancer_emergence_rate, :random_disruption_rate, :random_emergence_rate,
                  :disruption_significance_uncorrected, :emergence_significance_uncorrected,
                  :gene, :underfitted,
                  :disruption_significance, :emergence_significance,
                  :disruption_significance_uncorrected_fitting_aware, :emergence_significance_uncorrected_fitting_aware,
                  ) do

  def self.from_string(line)
      motif, official_gene_name, quality,
          cancer_to_random_disruption_ratio, disruption_significance_fitting_aware,
          cancer_to_random_emergence_ratio, emergence_significance_fitting_aware,
          random_disrupted, random_emerged, random_total_before_substitution, random_total_after_substitution,
          cancer_disrupted, cancer_emerged, cancer_total_before_substitution, cancer_total_after_substitution,
          cancer_disruption_rate, cancer_emergence_rate, random_disruption_rate, random_emergence_rate,
          disruption_significance_uncorrected, emergence_significance_uncorrected,
          gene, underfitted,
          disruption_significance, emergence_significance,
          disruption_significance_uncorrected_fitting_aware, emergence_significance_uncorrected_fitting_aware = line.chomp.split("\t")
      RateComparisonInfo.new(
          motif, official_gene_name, quality.upcase.to_sym,
          to_float(cancer_to_random_disruption_ratio), to_float(disruption_significance_fitting_aware),
          to_float(cancer_to_random_emergence_ratio), to_float(emergence_significance_fitting_aware),
          random_disrupted.to_i, random_emerged.to_i, random_total_before_substitution.to_i, random_total_after_substitution.to_i,
          cancer_disrupted.to_i, cancer_emerged.to_i, cancer_total_before_substitution.to_i, cancer_total_after_substitution.to_i,
          cancer_disruption_rate.to_f, cancer_emergence_rate.to_f, random_disruption_rate.to_f, random_emergence_rate.to_f,
          to_float(disruption_significance_uncorrected), to_float(emergence_significance_uncorrected),
          gene, underfitted,
          to_float(disruption_significance), to_float(emergence_significance),
          to_float(disruption_significance_uncorrected_fitting_aware), to_float(emergence_significance_uncorrected_fitting_aware),
        )
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.readline # skip header
      each_in_stream(f, &block)
    end
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.map{|line| RateComparisonInfo.from_string(line) }.each(&block)
  end

  def self.table_header
    RateComparisonInfo::COLUMN_ORDER.map{|column_id| RateComparisonInfo::COLUMN_NAMES[column_id] }.join("\t")
  end

  def to_s
    RateComparisonInfo::COLUMN_ORDER.map{|column_id| self[column_id] }.join("\t")
  end
end

RateComparisonInfo::COLUMN_NAMES = {
  motif: 'Motif',
  quality: 'Motif quality',
  gene: 'TF gene',
  official_gene_name: 'TF gene (official name)',
  random_disrupted: 'Random disrupted',
  random_emerged: 'Random emerged',
  random_total_before_substitution: 'Random total sites before subtitution',
  random_total_after_substitution: 'Random total sites after subtitution',
  cancer_disrupted: 'Cancer disrupted',
  cancer_emerged: 'Cancer emerged',
  cancer_total_before_substitution: 'Cancer total sites before subtitution',
  cancer_total_after_substitution: 'Cancer total sites after subtitution',
  cancer_disruption_rate: 'Cancer disruption rate',
  cancer_emergence_rate: 'Cancer emergence rate',
  random_disruption_rate: 'Random disruption rate',
  random_emergence_rate: 'Random emergence rate',
  cancer_to_random_disruption_ratio: 'Cancer to random disruption rate ratio',
  cancer_to_random_emergence_ratio: 'Cancer to random emergence rate ratio',
  disruption_significance_uncorrected: 'Significance of difference in disruption rate (not corrected)',
  emergence_significance_uncorrected: 'Significance of difference in emergence rate (not corrected)',
  disruption_significance: 'Significance of difference in disruption rate (corrected on multiple comparison)',
  emergence_significance: 'Significance of difference in emergence rate (corrected on multiple comparison)',
  underfitted: 'Number of sites underfitted',
  disruption_significance_uncorrected_fitting_aware: 'Significance of difference in disruption rate (not corrected, fitting aware)',
  emergence_significance_uncorrected_fitting_aware: 'Significance of difference in emergence rate (not corrected, fitting aware)',
  disruption_significance_fitting_aware: 'Significance of difference in disruption rate (corrected on multiple comparison, fitting aware)',
  emergence_significance_fitting_aware: 'Significance of difference in emergence rate (corrected on multiple comparison, fitting aware)',
}

RateComparisonInfo::COLUMN_ORDER = [:motif, :official_gene_name, :quality,
                  :cancer_to_random_disruption_ratio, :disruption_significance_fitting_aware,
                  :cancer_to_random_emergence_ratio, :emergence_significance_fitting_aware,
                  :random_disrupted, :random_emerged, :random_total_before_substitution, :random_total_after_substitution,
                  :cancer_disrupted, :cancer_emerged, :cancer_total_before_substitution, :cancer_total_after_substitution,
                  :cancer_disruption_rate, :cancer_emergence_rate, :random_disruption_rate, :random_emergence_rate,
                  :disruption_significance_uncorrected, :emergence_significance_uncorrected,
                  :gene, :underfitted,
                  :disruption_significance, :emergence_significance,
                  :disruption_significance_uncorrected_fitting_aware, :emergence_significance_uncorrected_fitting_aware,
                ]
