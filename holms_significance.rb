$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'


with_header = ARGV.delete('--with-header')
with_header_column = ARGV.delete('--with-header-column')
input_file = ARGV[0]

unless input_file
  raise "Specify input file with data in format:\n" +
        "disrupted_shuffle disrupted_cancer total_shuffle total_cancer\n" + 
        "..."
end


lines = File.readlines(input_file)
header = lines.shift  if with_header
data = lines.map{|line| line.strip.split("\t").map(&:to_i) }

if with_header_column
  header_column = data.map(&:first)
  counts = data.map{|line| line.drop(1) }
else
  counts = data
end


pvalues = data.map{|disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer|
  calculate_fisher_test(disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer)
}

holms_corrections = calculate_holms_correction(pvalues)

if with_header_column
  holms_corrections.zip(header_column).each do |pvalue_corrected, header_column_value|
    puts "#{header_column_value}\t#{pvalue_corrected}"
  end
else
  holms_corrections.each do |pvalue_corrected|
    puts pvalue_corrected
  end
end