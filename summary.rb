$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'
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

pvalue_correction_method = 'fdr'

OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <cancer statistics folder> <random stats folder> <motif names> <hocomoco gene infos> [options]"
  opts.separator('Options:')
  opts.on('--correction METHOD', 'P-value correction method (holm/fdr/hochberg/hommel/bonferroni/BH/BY/none -- it\'s processed by R). Default is fdr.') {|value|
    pvalue_correction_method = value
  }
end.parse!(ARGV)

raise 'Specify folder for cancer statistics'  unless cancer_dirname = ARGV[0] # './results/motif_statistics/cpg/cancer'
raise 'Specify folder for random statistics'  unless random_dirname = ARGV[1] # './results/motif_statistics/cpg/random'
raise 'Specify file with motif names'  unless motif_names_filename = ARGV[2] # './source_data/motif_names.txt'
raise 'Specify file with motif collection infos'  unless hocomoco_motifs_filename = ARGV[3] # './source_data/hocomoco_genes_infos.csv'

motif_names = File.readlines(motif_names_filename).map(&:strip)
motif_collection_infos = load_motif_infos(hocomoco_motifs_filename)

column_names =  {
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
}

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
  motif_infos[:cancer_disruption_rate][motif] = motif_infos[:cancer_disrupted][motif].to_f / motif_infos[:cancer_total_before_substitution][motif]
  motif_infos[:cancer_emergence_rate][motif] = motif_infos[:cancer_emerged][motif].to_f / motif_infos[:cancer_total_after_substitution][motif]

  motif_infos[:random_disruption_rate][motif] = motif_infos[:random_disrupted][motif].to_f / motif_infos[:random_total_before_substitution][motif]
  motif_infos[:random_emergence_rate][motif] = motif_infos[:random_emerged][motif].to_f / motif_infos[:random_total_after_substitution][motif]

  motif_infos[:cancer_to_random_disruption_ratio][motif] = motif_infos[:random_disruption_rate][motif] != 0 ? motif_infos[:cancer_disruption_rate][motif].to_f / motif_infos[:random_disruption_rate][motif] : nil
  motif_infos[:cancer_to_random_emergence_ratio][motif] = motif_infos[:random_emergence_rate][motif] != 0 ? motif_infos[:cancer_emergence_rate][motif].to_f / motif_infos[:random_emergence_rate][motif] : nil

  motif_infos[:disruption_significance_uncorrected][motif] = significance_calculator.calculate([motif_infos[:cancer_disrupted][motif], motif_infos[:random_disrupted][motif], motif_infos[:cancer_total_before_substitution][motif], motif_infos[:random_total_before_substitution][motif]])
  motif_infos[:emergence_significance_uncorrected][motif] = significance_calculator.calculate([motif_infos[:cancer_emerged][motif], motif_infos[:random_emerged][motif], motif_infos[:cancer_total_after_substitution][motif], motif_infos[:random_total_after_substitution][motif]])
}

motif_infos[:disruption_significance] = significance_corrector.correct_hash(motif_infos[:disruption_significance_uncorrected])
motif_infos[:emergence_significance] = significance_corrector.correct_hash(motif_infos[:emergence_significance_uncorrected])

motif_infos.merge!(motif_collection_infos)

output_columns = [:motif, :official_gene_name, :quality,
                  :cancer_to_random_disruption_ratio, :disruption_significance,
                  :cancer_to_random_emergence_ratio, :emergence_significance,
                  :random_disrupted, :random_emerged, :random_total_before_substitution, :random_total_after_substitution,
                  :cancer_disrupted, :cancer_emerged, :cancer_total_before_substitution, :cancer_total_after_substitution,
                  :cancer_disruption_rate, :cancer_emergence_rate, :random_disruption_rate, :random_emergence_rate,
                  :disruption_significance_uncorrected, :emergence_significance_uncorrected,
                  :gene]

puts output_columns.map{|column_id| column_names[column_id] }.join("\t")
motif_names.map{|motif_name|
  puts output_columns.map{|column_id| motif_infos[column_id][motif_name]}.join("\t")
}
