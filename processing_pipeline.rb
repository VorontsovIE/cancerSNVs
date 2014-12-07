$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'
require 'table'
require 'table_combiner'

def load_motif_qualities(filename)
  File.readlines(filename).drop(1).map{|line| motif, quality = line.chomp.split("\t").first(2); [motif.upcase, quality.upcase] }.to_h
end

def load_motif_infos(filename)
  table = Table.read(filename, with_header_row: true, with_header_column: true)
  table.header_row = { #name: 'NAME',   ## it stands in header
    gene: 'Gene',
    motif_quality: 'Motif quality',
    motif_weight: 'Motif weight',
    uniprot_human_tf: 'Human TFs',
    uniprot_mouse_tf: 'Mouse TFs',
    motif_consensus: 'Consensus',
  }
  table.select_columns([:gene, :motif_quality])
end

def load_count_table(filename)
  count_table = Table.read(filename, with_header_row: false, with_header_column: true)
  count_table.each{|line| line.map!(&:to_i) }
  count_table
end

def combine_counts(mutation_type, context_type, shuffle_type_filename, cancer_type_filename)
  files = ["results/disrupted/#{shuffle_type_filename}/#{mutation_type}_#{context_type}_disrupted_shuffle.txt",
          "results/disrupted/#{cancer_type_filename}/#{mutation_type}_#{context_type}_disrupted_cancer.txt",
          "results/disrupted/#{shuffle_type_filename}/#{mutation_type}_#{context_type}_total_shuffle.txt",
          "results/disrupted/#{cancer_type_filename}/#{mutation_type}_#{context_type}_total_cancer.txt"]

  count_tables = files.map{|filename| load_count_table(filename) }
  counts_table = TableColumnCombiner.new(count_tables).combine
  counts_table.header_row = Table::Header.create({shuffle_disrupted: 'Shuffle disrupted',
                                                  cancer_disrupted: 'Cancer disrupted',
                                                  shuffle_total: 'Shuffle total',
                                                  cancer_total: 'Cancer total'})

  pvalue_calculator = PvalueCalculator.new(class_counts: :class_and_total)
  pvalue_corrector = HolmsPvalueCorrector.new
  corrected_pvalues_table = calculate_corrected_pvalues(counts_table, pvalue_calculator, pvalue_corrector, header: [:significance, 'Corrected P-value'])
  
  TableColumnCombiner.new([counts_table, corrected_pvalues_table]).combine
end


def comparison_to_shuffle(mutation_type, context_type, motif_infos, shuffle_type_filename, cancer_type_filename)
  counts_table_w_pvalues = combine_counts(mutation_type, context_type, shuffle_type_filename, cancer_type_filename)
  File.open("results/disrupted/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv", 'w') do |fw|
    counts_table_w_pvalues.output(fw)
  end

  counts_table_w_pvalues.add_columns({ disruption_rate_shuffle: 'Disruption rate (shuffle)',
                                      disruption_rate_cancer: 'Disruption rate (cancer)',
                                      cancer_to_shuffle_rate: 'Cancer to shuffle disruption rate',
                                    })
   counts_table_w_pvalues.each_row do |row|
    motif = row.row_name
    row[:disruption_rate_shuffle] = row[:shuffle_disrupted].to_f / row[:shuffle_total]
    row[:disruption_rate_cancer] = row[:cancer_disrupted].to_f / row[:cancer_total]
    row[:cancer_to_shuffle_rate] = row[:disruption_rate_cancer] / row[:disruption_rate_shuffle]
  end

  table = TableColumnCombiner.new([counts_table_w_pvalues, motif_infos]).combine
  table.add_column({official_gene_name: 'Official gene name'})
  table.each_row do |row|
    row[:official_gene_name] = row[:gene].split.first
    row[:gene] = row[:gene].split.join(',')
  end
  File.open("results/disrupted/augmented/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv", 'w') do |fw|
    table.output(fw)
  end

  table_selected = table.select_rows{|row|
    row[:significance] <= 0.05 && row[:motif_quality] != 'D'
  }.sort_rows_by{|row|
    - row[:cancer_to_shuffle_rate]
  }
  File.open("results/disrupted/filtered/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv", 'w') do |fw|
    table_selected.output(fw)
  end
end

['results/disrupted/', 'results/disrupted/filtered', 'results/disrupted/augmented'].each{|dirname| Dir.mkdir(dirname)  unless Dir.exist?(dirname) }

# system "ruby calculate_motif_statistics.rb source_data/cancer_SNPs.txt results/disrupted/cancer/cancer.txt"
# system "ruby calculate_motif_statistics.rb source_data/shuffle_SNPs.txt results/disrupted/shuffle/shuffle.txt"

motif_infos = load_motif_infos('source_data/hocomoco_genes_infos.csv')

# ['shuffle', 'sampled_shuffle_log2', 'sampled_shuffle_log10', 'sampled_shuffle_0.5log10', 'contexted_sampled_shuffle_log10'].each do |shuffle_type_filename|
#   # [:regulatory, :intronic, :promoter].each do |mutation_type|
#   [:regulatory].each do |mutation_type|
#     [:cpg, :tpc, :not_cpg_tpc, :any_context].each do |context_type|
#       comparison_to_shuffle(mutation_type, context_type, motif_infos, shuffle_type_filename, 'cancer')
#     end
#   end
# end

[:regulatory].each do |mutation_type|
  [:cpg, :tpc, :not_cpg_tpc, :any_context].each do |context_type|
    comparison_to_shuffle(mutation_type, context_type, motif_infos, 'shuffle_created', 'cancer_created')
  end
end
