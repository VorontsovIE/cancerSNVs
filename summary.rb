$:.unshift File.absolute_path('lib', __dir__)
require 'statistics/statistical_significance'
require 'rate_comparison_infos'
require 'statistics/fisher_table'
require 'fitting/fitting_logs'
require 'uniprot_info'
require 'optparse'
require_relative 'experiment_configuration'

def read_motif_counts(filename)
  motif_counts = File.readlines(filename).map{|line|
    motif_name, count = line.chomp.split("\t")
    [motif_name, count.to_i]
  }.to_h
end

def load_motif_qualities(filename)
  File.readlines(filename).map{|line|
    motif, quality = line.chomp.split("\t")
    [motif, quality.to_sym]
  }.to_h
end

motif_qualities_to_take = [:A, :B, :C]
pvalue_correction_method = 'fdr'
control_set_multiplier = 1
ignore_underfitting = false
fitting_log_filename = nil

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <cancer statistics folder> <random stats folder> --fitting-log <fitting_log> [options]"
  opts.separator('Options:')
  opts.on('--correction METHOD', 'P-value correction method (holm/fdr/hochberg/hommel/bonferroni/BH/BY/none -- it\'s processed by R). Default is fdr.') {|value|
    pvalue_correction_method = value
  }
  opts.on('--expand-control-set N', 'Calculate statistics as if numbers in control set were N times greater. Use it only for preliminary checks') {|value|
    control_set_multiplier = Integer(value)
  }
  opts.on('--ignore-underfitting', 'Don\'t take underfitted values into account. Use it only for preliminary checks' ) {
    ignore_underfitting = true
  }
  opts.on('--fitting-log FILE', 'Specify fitting log file to make underfitting corrections') {|filename|
    fitting_log_filename = filename
  }
  opts.on('--motif-qualities QUALITIES', 'Take motifs of given qualities (comma-separated). Default is A,B,C') {|value|
    motif_qualities_to_take = value.upcase.split(',').map(&:to_sym)
  }
end.parse!(ARGV)

raise 'Specify folder for cancer statistics'  unless cancer_dirname = ARGV[0] # './results/motif_statistics/cpg/cancer'
raise 'Specify folder for random statistics'  unless random_dirname = ARGV[1] # './results/motif_statistics/cpg/random'

motif_names = File.readlines(LocalPaths::Secondary::MotifNames).map(&:strip)
motif_qualities = load_motif_qualities(LocalPaths::Secondary::MotifQualities)

uniprot_infos_by_motif = read_uniprot_infos_by_motif(LocalPaths::Secondary::HocomocoUniprots, LocalPaths::Secondary::UniprotDump)

fitting_logs = load_motif_underfitting_rates(fitting_log_filename)

motif_infos = {
  random_disrupted: File.join(random_dirname, "sites_disrupted.txt"),
  random_emerged: File.join(random_dirname, "sites_emerged.txt"),
  random_total_before_substitution: File.join(random_dirname, "sites_before.txt"),
  random_total_after_substitution: File.join(random_dirname, "sites_after.txt"),
  random_core: File.join(random_dirname, "substitutions_in_core.txt"),
  random_flank: File.join(random_dirname, "substitutions_in_flank.txt"),

  cancer_disrupted: File.join(cancer_dirname, "sites_disrupted.txt"),
  cancer_emerged: File.join(cancer_dirname, "sites_emerged.txt"),
  cancer_total_before_substitution: File.join(cancer_dirname, "sites_before.txt"),
  cancer_total_after_substitution: File.join(cancer_dirname, "sites_after.txt"),
  cancer_core: File.join(cancer_dirname, "substitutions_in_core.txt"),
  cancer_flank: File.join(cancer_dirname, "substitutions_in_flank.txt"),
}.map {|column_name, filename|
  [column_name, read_motif_counts(filename)]
}.to_h

significance_corrector = PvalueCorrector.new(pvalue_correction_method)

motif_statistics = motif_names.select{|motif|
  motif_qualities_to_take.include?( motif_qualities[motif] )
}.map{|motif|
  MotifStatistics.new(
    motif: motif,
    disruption_table: FisherTable.by_class_and_total(
      class_a_total: motif_infos[:cancer_total_before_substitution][motif],
      class_a_positive: motif_infos[:cancer_disrupted][motif],
      class_b_total: motif_infos[:random_total_before_substitution][motif] * control_set_multiplier,
      class_b_positive: motif_infos[:random_disrupted][motif] * control_set_multiplier
    ),
    emergence_table: FisherTable.by_class_and_total(
      class_a_total: motif_infos[:cancer_total_after_substitution][motif],
      class_a_positive: motif_infos[:cancer_emerged][motif],
      class_b_total: motif_infos[:random_total_after_substitution][motif] * control_set_multiplier,
      class_b_positive: motif_infos[:random_emerged][motif] * control_set_multiplier
    ),

    core_flank_table: FisherTable.by_two_classes(
      class_a_positive: motif_infos[:cancer_core][motif],
      class_a_negative: motif_infos[:cancer_flank][motif],
      class_b_positive: motif_infos[:random_core][motif],
      class_b_negative: motif_infos[:random_flank][motif]
    ),
    random_unclassified: ignore_underfitting ? 0 : fitting_logs.fetch(motif, 0) * control_set_multiplier,

    gene: uniprot_infos_by_motif[motif].map(&:all_gene_names).uniq.join(','),
    quality: motif_qualities[motif],
    official_gene_name: uniprot_infos_by_motif[motif].map(&:primary_gene_name).uniq.join(',')
  )
}

puts MotifCollectionStatistics.new(motif_statistics, pvalue_corrector: significance_corrector).to_s
