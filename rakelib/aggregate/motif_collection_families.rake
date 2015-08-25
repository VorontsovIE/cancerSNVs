desc 'Extract motif collection table'
task :motif_collection_families do
	motifs = File.readlines(LocalPaths::Secondary::MotifNames).map(&:strip)
	motif_qualities = load_motif_qualities(LocalPaths::Secondary::MotifQualities)

	File.open('results/motif_collection.tsv', 'w') do |fw|
		motifs.each do |motif|
			families_3 = MOTIF_FAMILY_RECOGNIZERS[3].subfamilies_by_motif(motif).map(&:to_s).join('; ')
			families_4 = MOTIF_FAMILY_RECOGNIZERS[4].subfamilies_by_motif(motif).map(&:to_s).join('; ')
			fw.puts [motif, motif_qualities[motif], families_3, families_4].join("\t")
		end
	end
end
