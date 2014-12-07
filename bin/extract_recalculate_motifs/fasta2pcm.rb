require 'shellwords'
LETTER_INDEX = {'A' => 0, 'C' =>1, 'G' => 2, 'T' => 3}

def pcm_from_fasta(seqs)
  seqs = seqs.map{|l| l.chomp.upcase }

  len = seqs.first.length

  pcm = Array.new(len){[0,0,0,0]}
  seqs.each do |seq|
    seq.each_char.each_with_index do |let,ind|
      pcm[ind][LETTER_INDEX[let]] += 1
    end
  end
  pcm
end

Shellwords.split($stdin.read).each do |filename|
  dirname = File.join(File.dirname(filename),'pcm')
  Dir.mkdir(dirname)  unless Dir.exist?(dirname)
  pcm = pcm_from_fasta(File.readlines(filename))
  File.open(File.join(dirname, File.basename(filename, File.extname(filename)) + '.pcm'), 'w') do |fw|
    fw.puts( pcm.map{|pos| pos.join("\t") }.join("\n") )
  end
end
