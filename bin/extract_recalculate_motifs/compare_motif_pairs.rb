motifs = File.readlines('./source_data/motif_names.txt').map(&:chomp)
dir_1 = './results/motif_sites/pcm'
dir_2 = './results/motif_sites_fitted/pcm'

motifs.each do |motif|
  mot_1 = File.join(dir_1, motif + '.pcm')
  mot_2 = File.join(dir_2, motif + '.pcm')
  result = `java -cp ape.jar ru.autosome.macroape.EvalSimilarity #{mot_1} #{mot_2} --pcm`.lines.map(&:chomp)
  sim = result.detect{|l| l.start_with?("S\t") }.split("\t").last
  puts "#{motif}\t#{sim}"
end
