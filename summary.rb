$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'
require 'rate_comparison_infos'
require 'optparse'

def read_motif_counts(filename)
  motif_counts = File.readlines(filename).map{|line|
    motif_name, count = line.chomp.split("\t")
    [motif_name, count.to_i]
  }.to_h
end

def load_motif_infos(filename)
  results = Hash.new{|h,k| h[k] = {} }
  File.readlines(filename).drop(1).each{|line|
    motif, gene, quality, weight, human_uniprot, mouse_uniprot, consensus = line.chomp.split("\t")
    results[:gene][motif] = gene.split.join(',')
    results[:quality][motif] = quality.upcase.to_sym
    results[:official_gene_name][motif] = gene.split.first
  }
  results
end

def load_motif_underfitting_rates(fitting_log_filename)
  return {}  unless File.exist?(fitting_log_filename)
  File.readlines(fitting_log_filename).drop(4).slice_before{|line|
    line.match(/^\t[^\t]/) && ! line.match(/(\d+) underfitted/)
  }.map{|lines|
    motif = lines[0].strip
    underfitted = lines[1].strip.match(/^(\d+) underfitted/)[1].to_i
    [motif, underfitted]
  }.to_h
end

pvalue_correction_method = 'fdr'
control_set_multiplier = 1
ignore_underfitting = false

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <cancer statistics folder> <random stats folder> <motif names> <hocomoco gene infos> [options]"
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
end.parse!(ARGV)

raise 'Specify folder for cancer statistics'  unless cancer_dirname = ARGV[0] # './results/motif_statistics/cpg/cancer'
raise 'Specify folder for random statistics'  unless random_dirname = ARGV[1] # './results/motif_statistics/cpg/random'
raise 'Specify file with motif names'  unless motif_names_filename = ARGV[2] # './source_data/motif_names.txt'
raise 'Specify file with motif collection infos'  unless hocomoco_motifs_filename = ARGV[3] # './source_data/hocomoco_genes_infos.csv'
raise 'Specify fitting log file'  unless fitting_log_filename = ARGV[4] # './results/motif_statistics/fitting_log/any/random_genome_13.log'

motif_names = File.readlines(motif_names_filename).map(&:strip)
motif_collection_infos = load_motif_infos(hocomoco_motifs_filename)
fitting_logs = load_motif_underfitting_rates(fitting_log_filename)

column_names =  RateComparisonInfo::COLUMN_NAMES

motif_infos = {
  random_disrupted: File.join(random_dirname, "sites_disrupted.txt"),
  random_emerged: File.join(random_dirname, "sites_emerged.txt"),
  random_total_before_substitution: File.join(random_dirname, "sites_before.txt"),
  random_total_after_substitution: File.join(random_dirname, "sites_after.txt"),

  cancer_disrupted: File.join(cancer_dirname, "sites_disrupted.txt"),
  cancer_emerged: File.join(cancer_dirname, "sites_emerged.txt"),
  cancer_total_before_substitution: File.join(cancer_dirname, "sites_before.txt"),
  cancer_total_after_substitution: File.join(cancer_dirname, "sites_after.txt"),
}.map {|column_name, filename|
  [column_name, read_motif_counts(filename)]
}.to_h

motif_infos.default_proc = ->(hsh,k) { hsh[k] = {} }

significance_calculator = PvalueCalculator.new(class_counts: :class_and_total)
significance_corrector = PvalueCorrector.new(pvalue_correction_method)

motif_names.each {|motif|
  motif_infos[:motif][motif] = motif

  cancer_disrupted = motif_infos[:cancer_disrupted][motif]
  cancer_total_before = motif_infos[:cancer_total_before_substitution][motif]
  cancer_emerged = motif_infos[:cancer_emerged][motif]
  cancer_total_after = motif_infos[:cancer_total_after_substitution][motif]

  random_disrupted = motif_infos[:random_disrupted][motif] * control_set_multiplier
  random_total_before = motif_infos[:random_total_before_substitution][motif] * control_set_multiplier
  random_emerged = motif_infos[:random_emerged][motif] * control_set_multiplier
  random_total_after = motif_infos[:random_total_after_substitution][motif] * control_set_multiplier

  cancer_disruption_rate = cancer_disrupted.to_f / cancer_total_before
  random_disruption_rate = random_disrupted.to_f / random_total_before

  cancer_emergence_rate = cancer_emerged.to_f / cancer_total_after
  random_emergence_rate = random_emerged.to_f / random_total_after

  random_underfitted = ignore_underfitting ? 0 : fitting_logs.fetch(motif, 0) * control_set_multiplier

  ###
  motif_infos[:cancer_disruption_rate][motif] = cancer_disruption_rate
  motif_infos[:cancer_emergence_rate][motif] = cancer_emergence_rate

  motif_infos[:random_disruption_rate][motif] = random_disruption_rate
  motif_infos[:random_emergence_rate][motif] = random_emergence_rate

  motif_infos[:cancer_to_random_disruption_ratio][motif] = random_disruption_rate != 0 ? cancer_disruption_rate / random_disruption_rate : nil
  motif_infos[:cancer_to_random_emergence_ratio][motif] = random_emergence_rate != 0 ? cancer_emergence_rate / random_emergence_rate : nil

  motif_infos[:disruption_significance_uncorrected][motif] = significance_calculator.calculate([cancer_disrupted, random_disrupted, cancer_total_before, random_total_before])
  motif_infos[:emergence_significance_uncorrected][motif] = significance_calculator.calculate([cancer_emerged, random_emerged, cancer_total_after, random_total_after])

  motif_infos[:underfitted][motif] = random_underfitted

  motif_infos[:disruption_significance_uncorrected_fitting_aware][motif] = significance_calculator.calculate([cancer_disrupted, random_disrupted, cancer_total_before, random_total_before], unclassified_a: 0, unclassified_b: random_underfitted)
  motif_infos[:emergence_significance_uncorrected_fitting_aware][motif] = significance_calculator.calculate([cancer_emerged, random_emerged, cancer_total_after, random_total_after], unclassified_a: 0, unclassified_b: random_underfitted)
}

motif_infos[:disruption_significance] = significance_corrector.correct_hash(motif_infos[:disruption_significance_uncorrected])
motif_infos[:emergence_significance] = significance_corrector.correct_hash(motif_infos[:emergence_significance_uncorrected])

motif_infos[:disruption_significance_fitting_aware] = significance_corrector.correct_hash(motif_infos[:disruption_significance_uncorrected_fitting_aware])
motif_infos[:emergence_significance_fitting_aware] = significance_corrector.correct_hash(motif_infos[:emergence_significance_uncorrected_fitting_aware])

motif_infos.merge!(motif_collection_infos)

output_columns = RateComparisonInfo::COLUMN_ORDER

puts RateComparisonInfo.table_header
motif_names.map{|motif_name|
  puts output_columns.map{|column_id| motif_infos[column_id][motif_name]}.join("\t")
}
