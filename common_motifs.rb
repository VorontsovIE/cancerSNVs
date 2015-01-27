require 'set'
require 'shellwords'

if $stdin.tty?
  files = ARGV
else
  files = Shellwords.split($stdin.read)
end

motif_sets = files.map do |fn|
  File.readlines(fn).drop(1).map{|l| l.chomp.split("\t").first }.to_set
end

puts motif_sets.inject(&:intersection).sort.join(', ')
