$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'
require 'table'
require 'optparse'

pvalue_calculator = PvalueCalculator.new(class_counts: :class_and_total)
pvalue_corrector = HolmsPvalueCorrector.new

with_header_row = false
with_header_column = false

OptionParser.new do |opts|
  opts.on('--no-correction', 'Don\'t correct pvalues') { pvalue_corrector = IdlePvalueCorrector.new }

  opts.on('--class-and-total-counts', 'Counts should be either in two_classes or class_and_total (default) format') { |value|
    pvalue_calculator = PvalueCalculator.new(class_counts: :class_and_total)
  }
  opts.on('--two-classes-counts', 'Counts should be either in two_classes or class_and_total (default) format') { |value|
    pvalue_calculator = PvalueCalculator.new(class_counts: :two_classes)
  }

  opts.on('--with-header', 'Input file has header row') { with_header_row = true }
  opts.on('--with-header-column', 'Input file has header column') { with_header_column = true }
end.parse!(ARGV)

table = Table.read(ARGF, with_header_row: with_header_row, with_header_column: with_header_column)
pvalues = table.map do |line|
  class_counts = line.map(&:to_i)
  raise "Bad input, should be 4 class count columns. See line:\n#{line.inspect}"  unless class_counts.size == 4
  pvalue_calculator.calculate(class_counts)
end

pvalues_corrected = pvalue_corrector.correct(pvalues)

table.add_column('P-value', pvalues_corrected)
table.output($stdout)
