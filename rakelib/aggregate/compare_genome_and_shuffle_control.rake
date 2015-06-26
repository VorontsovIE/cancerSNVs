directory 'results/motif_statistics/shuffle_vs_genome_control'
desc 'Collect sample statistics'
task :compare_genome_and_shuffle_control => 'results/motif_statistics/shuffle_vs_genome_control' do
  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      results = []
      Configuration.getAlexandrovWholeGenomeCancers.each do |cancer_type|
        folder = File.join('results/motif_statistics/common/', 'Alexandrov', cancer_type.to_s, 'any', protected_or_subjected.to_s, characteristic.to_s)
        prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
        motifs_determined_by_shuffle_control = File.readlines(File.join(folder, 'compared_to_each_shuffle.txt'))
        motifs_determined_by_genome_control  = File.readlines(File.join(folder, 'compared_to_each_genome.txt'))
        motifs_determined_by_both_controls   = File.readlines(File.join(folder, 'compared_to_each.txt'))

        results << [cancer, motifs_determined_by_shuffle_control.size, motifs_determined_by_genome_control.size, motifs_determined_by_both_controls.size]
      end
      output_filename = File.join('results/motif_statistics/shuffle_vs_genome_control', "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.tsv")
      File.write(output_filename, results.map{|row| row.join("\t") }.join("\n"))
    end
  end
end

task :default => [:compare_genome_and_shuffle_control]
