$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'

['results/disrupted/', 'results/disrupted/filtered', 'results/disrupted/augmented'].each{|dirname| Dir.mkdir(dirname)  unless Dir.exist?(dirname) }

# system "ruby calculate_motif_statistics.rb source_data/cancer_SNPs.txt results/disrupted/cancer/cancer.txt"
# system "ruby calculate_motif_statistics.rb source_data/shuffle_SNPs.txt results/disrupted/shuffle/shuffle.txt"

def load_motif_qualities(filename)
  File.readlines(filename).drop(1).map{|line| motif, quality = line.chomp.split("\t").first(2); [motif.upcase, quality.upcase] }.to_h
end

def combine_results(mutation_type, context_type, shuffle_type_filename)
  IO.popen(["ruby",  "combine_files.rb",
                    "results/disrupted/#{shuffle_type_filename}/#{mutation_type}_#{context_type}_disrupted_shuffle.txt",
                    "results/disrupted/cancer/#{mutation_type}_#{context_type}_disrupted_cancer.txt",
                    "results/disrupted/#{shuffle_type_filename}/#{mutation_type}_#{context_type}_total_shuffle.txt",
                    "results/disrupted/cancer/#{mutation_type}_#{context_type}_total_cancer.txt",
                    "--header", "one"]) do |counts_io|
    with_temp_file("counts.txt") do |counts_fw|
      counts_fw.puts("motif\tdisrupted shuffle\tdisrupted cancer\ttotal shuffle\ttotal cancer")
      counts_fw.print(counts_io.read)
      counts_fw.close

      IO.popen("ruby holms_significance.rb --with-header --with-header-column #{counts_fw.path}") do |significances_io|
        with_temp_file("significances.txt") do |significances_fw|
          significances_fw.puts("motif\tHolms corrected P-value")
          significances_fw.print(significances_io.read)
          significances_fw.close

          IO.popen(["ruby", "combine_files.rb", counts_fw.path, significances_fw.path, "--header", "one"]) do |result_io|
            File.write("results/disrupted/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv", result_io.read)
          end
        end
      end
    end
  end
end

genes_by_motif = File.readlines('source_data/hocomoco_genes.csv').drop(1)
                      .map{|line| line.chomp.split("\t").first(2) }
                      .map{|motif,genes| [motif.sub(/\bHOCOMOCO\/(.+)\.pat\b/, '\1'), genes.split.map(&:strip)] }.to_h


motif_qualities = load_motif_qualities('source_data/hocomoco_info_AD.csv')

['shuffle', 'sampled_shuffle_log2', 'sampled_shuffle_log10', 'sampled_shuffle_0.5log10', 'contexted_sampled_shuffle_log10'].each do |shuffle_type_filename|
# [:regulatory, :intronic, :promoter].each do |mutation_type|
[:regulatory].each do |mutation_type|
  [:cpg, :tpc, :not_cpg_tpc, :any_context].each do |context_type|
    combine_results(mutation_type, context_type, shuffle_type_filename)
    lines = File.readlines("results/disrupted/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv")
    infos = lines.drop(1).map do |line|
      motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance = line.chomp.split("\t")
      disrupted_shuffle = disrupted_shuffle.to_i
      disrupted_cancer = disrupted_cancer.to_i
      total_shuffle = total_shuffle.to_i
      total_cancer = total_cancer.to_i
      holms_significance = holms_significance.to_f
      disruption_rate_shuffle = disrupted_shuffle.to_f / total_shuffle.to_f
      disruption_rate_cancer = disrupted_cancer.to_f / total_cancer.to_f
      cancer_to_shuffle_rate = disruption_rate_cancer / disruption_rate_shuffle
      genes = genes_by_motif[motif.upcase.tr('+!','__')].join(',')
      genes_official = genes_by_motif[motif.upcase.tr('+!','__')].first
      [ motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_qualities[motif.upcase]]
    end

    File.open("results/disrupted/augmented/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv", 'w') do |fw|
      fw.puts [lines.first.chomp, 'shuffle disruption rate', 'cancer disruption rate', 'cancer to shuffle rate', 'genes', 'gene official symbol', 'motif_quality'].join("\t")
      infos.sort_by do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality|
        motif
      end.each do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality|
        fw.puts [motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality].join("\t")
      end
    end

    File.open("results/disrupted/filtered/#{mutation_type}_#{context_type}_#{shuffle_type_filename}.csv", 'w') do |fw|
      fw.puts [lines.first.chomp, 'shuffle disruption rate', 'cancer disruption rate', 'cancer to shuffle rate', 'genes', 'gene official symbol', 'motif_quality'].join("\t")
      infos.select do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality|
        holms_significance <= 0.05 && motif_qualities[motif.upcase] != 'D' # && motif_qualities[motif.upcase] != 'C'
      end.sort_by do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality|
        cancer_to_shuffle_rate
      end.each do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality|
        fw.puts [motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes, genes_official, motif_quality].join("\t")
      end
    end

  end
end
end