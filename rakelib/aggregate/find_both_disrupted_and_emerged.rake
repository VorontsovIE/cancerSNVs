require 'set'
require_relative '../../lib/motif_family_recognizer'


directory 'results/motif_statistics/disruption_and_emergence'
desc 'Find motifs subjected or protected from both disruption and emergence'
task :find_both_disrupted_and_emerged => 'results/motif_statistics/disruption_and_emergence' do
  motif_qualities = load_motif_qualities(LocalPaths::Secondary::GeneInfos)

  [:protected, :subjected].each do |protected_or_subjected|
    prep = (protected_or_subjected == :subjected) ? 'to' : 'from'

    both_disrupted_and_emerged_by_sample = Configuration.sample_paths(with_nik_zainal: false, with_yeast: false).map{|sample_name, sample_folder|
      disruption_motifs_fn = File.join('results/motif_statistics/common', sample_folder, 'any', protected_or_subjected.to_s, 'disruption/compared_to_each.txt')
      emergence_motifs_fn = File.join('results/motif_statistics/common', sample_folder, 'any', protected_or_subjected.to_s, 'emergence/compared_to_each.txt')
      disruption_motifs = File.readlines(disruption_motifs_fn).map(&:strip)
      emergence_motifs = File.readlines(emergence_motifs_fn).map(&:strip)
      both_disrupted_and_emerged = Set.new(disruption_motifs) & Set.new(emergence_motifs)
      [sample_name, both_disrupted_and_emerged]
    }.to_h

    all_motifs = both_disrupted_and_emerged_by_sample.map{|sample_name, motifs| motifs }.inject(Set.new, :|).sort
    all_samples = Configuration.sample_paths(with_nik_zainal: false, with_yeast: false).keys

    results = []
    results << ['Motif', 'Motif quality', 'Motif families (level 3)', 'Motif families (level 4)', *all_samples]
    all_motifs.each{|motif|
      quality = motif_qualities[motif]
      families_3 = MOTIF_FAMILY_RECOGNIZERS[3].subfamilies_by_motif(motif).map(&:to_s).join('; ')
      families_4 = MOTIF_FAMILY_RECOGNIZERS[4].subfamilies_by_motif(motif).map(&:to_s).join('; ')
      occurences_in_sample = all_samples.map{|sample| both_disrupted_and_emerged_by_sample[sample].include?(motif) ? 1 : nil }
      results << [motif, quality, families_3, families_4, *occurences_in_sample]
    }

    output_filename = File.join('results/motif_statistics/disruption_and_emergence', "#{protected_or_subjected}_#{prep}_both_in_any_context.tsv")
    File.write(output_filename, results.map{|row| row.join("\t") }.join("\n"))
  end
end

task :default => [:find_both_disrupted_and_emerged]
