$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'

# system "ruby calculate_motif_statistics.rb source_data/cancer_SNPs.txt results/disrupted/cancer/cancer.txt"
# system "ruby calculate_motif_statistics.rb source_data/shuffle_SNPs.txt results/disrupted/shuffle/shuffle.txt"

def load_motif_qualities(filename)
  File.readlines(filename).drop(1).map{|line| motif, quality = line.chomp.split("\t").first(2); [motif.upcase, quality.upcase] }.to_h
end

def combine_results(mutation_type, context_type)
    IO.popen(["ruby",  "combine_files.rb",
                      "results/disrupted/shuffle/#{mutation_type}_#{context_type}_disrupted_shuffle.txt",
                      "results/disrupted/cancer/#{mutation_type}_#{context_type}_disrupted_cancer.txt",
                      "results/disrupted/shuffle/#{mutation_type}_#{context_type}_total_shuffle.txt",
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
              File.write("results/#{mutation_type}_#{context_type}.csv", result_io.read)
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

# [:regulatory, :intronic, :promoter].each do |mutation_type|
[:regulatory].each do |mutation_type|
  [:cpg, :tpc, :not_cpg_tpc, :any_context].each do |context_type|
    combine_results(mutation_type, context_type)
    lines = File.readlines("results/#{mutation_type}_#{context_type}.csv")
    File.open("results/filtered/#{mutation_type}_#{context_type}.csv", 'w') do |fw|
      fw.puts [lines.first.chomp, 'shuffle disruption rate', 'cancer disruption rate', 'cancer to shuffle rate', 'genes'].join("\t")
      lines.drop(1).map do |line|
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
        [ motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes ]
      end.select do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes|
        holms_significance <= 0.05 && motif_qualities[motif.upcase] != 'D' # && motif_qualities[motif.upcase] != 'C'
      end.sort_by do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes|
        cancer_to_shuffle_rate
      end.each do |motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes|
        fw.puts [motif, disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer, holms_significance, disruption_rate_shuffle, disruption_rate_cancer, cancer_to_shuffle_rate, genes].join("\t")
      end
    end
  end
end
