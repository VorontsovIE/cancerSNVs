# Convert Alexandrov et al. supplementary table 4 from xls into csv
require 'roo'
require 'roo-xls'

raise 'Specify xls file from Alexandrov et al. supplementary table 4'  unless supplementary_table_filename = ARGV[0]
xls = Roo::Excel.new(supplementary_table_filename)
xls.drop(3).map{|line|
  line.map{|el| el.is_a?(Numeric) ? el.to_i.to_s : el }.join("\t")
}.each{|x| puts x}
