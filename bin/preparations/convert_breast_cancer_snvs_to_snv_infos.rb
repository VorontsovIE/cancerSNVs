$:.unshift File.absolute_path('./../../lib', __dir__)
require 'data_import/breast_cancer_snv'
require 'snv_info'
require 'experiment_configuration'

raise 'Specify filename'  unless filename = ARGV.first

# It will output all mutations in SNVInfo format in pyrimidine context.
# Strand indicates strand of pyrimidine context
puts SNVInfo::HEADER
BreastCancerSNV.each_in_file(filename).map{|snv| 
  snv.to_snv_info(GENOME_READER, flank_length: 50)
}.map(&:in_pyrimidine_context)
.each{|snv_info| 
  puts snv_info
}
