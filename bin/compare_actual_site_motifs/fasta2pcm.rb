require 'shellwords'
require 'optparse'
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

def matrix_to_s(matrix)
  matrix.map{|pos| pos.join("\t") }.join("\n")
end

stdin_mode = :filelist

OptionParser.new do |opts|
  opts.on('--stdin-fasta', 'Use stdin as a source of FASTA. Output to stdout.') { stdin_mode = :fasta}
  opts.on('--stdin-filelist', 'Use stdin as a list of FASTA files. Output to a list of new files. This is default') { stdin_mode = :filelist }
end.parse!(ARGV)

if stdin_mode == :filelist
  Shellwords.split($stdin.read).each do |filename|
    dirname = File.join(File.dirname(filename),'pcm')
    Dir.mkdir(dirname)  unless Dir.exist?(dirname)
    pcm = pcm_from_fasta(File.readlines(filename))
    File.open(File.join(dirname, File.basename(filename, File.extname(filename)) + '.pcm'), 'w') do |fw|
      fw.puts(matrix_to_s(pcm))
    end
  end
elsif stdin_mode == :fasta
  pcm = pcm_from_fasta($stdin.readlines)
  puts(matrix_to_s(pcm))
else
  raise LogicError, "Shouldn't be here"
end
