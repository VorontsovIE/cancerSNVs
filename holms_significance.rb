$:.unshift File.absolute('lib', __dir__)
require 'statistical_significance'


with_header = ARGV.delete('--with-header')
input_file = ARGV[0]

unless input_file
  raise "Specify input file with data in format:\n" +
        "disrupted_shuffle disrupted_cancer total_shuffle total_cancer\n" + 
        "..."
end


lines = File.readlines(input_file)
lines.shift  if with_header
data = lines.map{|line| line.strip.split.map(&:to_i) }


pvalues = data.map{|disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer|
  calculate_fisher_test(disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer)
}

holms_corrections = calculate_holms_correction(pvalues)

holms_corrections.each do |h|
  puts h
end
