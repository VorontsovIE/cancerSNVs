require 'set'

dir_1 = ARGV[0] # './results/motif_sites/pcm'
dir_2 = ARGV[1] # './results/motif_sites_fitted/pcm'

motifs_names_1 = Dir.glob(File.join(dir_1,'*.pcm')).map{|fn| File.basename(fn, '.pcm') }.to_set
motifs_names_2 = Dir.glob(File.join(dir_2,'*.pcm')).map{|fn| File.basename(fn, '.pcm') }.to_set
motif_names = (motifs_names_1 & motifs_names_2).sort

motif_names.each do |motif|
  mot_1 = File.join(dir_1, motif + '.pcm')
  mot_2 = File.join(dir_2, motif + '.pcm')
  result = `java -cp ape.jar ru.autosome.macroape.EvalSimilarity #{mot_1} #{mot_2} --pcm`.lines.map(&:chomp)
  sim = result.detect{|l| l.start_with?("S\t") }.split("\t").last
  puts "#{motif}\t#{sim}"
end
