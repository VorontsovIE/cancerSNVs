$:.unshift File.absolute_path('../../lib', __dir__)
require 'snv_info'

raise 'Specify file with sequences'  unless sequences_filename = ARGV[0] # './source_data/sequences_background_genome.txt'

# Can be easily replaced with SNVInfo, because this script is not needed in this form
puts SNVInfo::HEADER
File.readlines(sequences_filename).each.lazy.map {|line|
  full_name, seq = line.chomp.split("\t")
  full_pos, mutated_to = full_name.split("/")
  chromosome, position = full_pos.split(':')
  position = position.to_i
  chromosome = chromosome.sub(/\Achr/, '')

  seq_w_snv = SequenceWithSNV.from_string(seq)
  strand = seq_w_snv.in_pyrimidine_context?  ?  :+  :  :-

  variant_id = full_pos
  sample_id = 'random genome mutations'
  SNVInfo.new(variant_id, seq_w_snv.in_pyrimidine_context,
              sample_id, chromosome, position, strand,
              RegionType.new)
}.each {|snv_info|
  puts snv_info.to_s
}
