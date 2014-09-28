$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'

# system "ruby calculate_motif_statistics.rb source_data/cancer_SNPs.txt results/disrupted/cancer/cancer.txt"
# system "ruby calculate_motif_statistics.rb source_data/shuffle_SNPs.txt results/disrupted/shuffle/shuffle.txt"

[:regulatory, :intronic, :promoter].each do |mutation_type|
  [:cpg, :tpc, :not_cpg_tpc, :any_context].each do |context_type|

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
              File.write("results/#{mutation_type}_#{context_type}.txt", result_io.read)
            end
          end
        end
      end
    end

  end
end

