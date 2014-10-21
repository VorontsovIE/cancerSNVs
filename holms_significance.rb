$:.unshift File.absolute_path('lib', __dir__)
require 'statistical_significance'
require 'table'
require 'table_combiner'
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

counts_table = Table.read(ARGF, with_header_row: with_header_row, with_header_column: with_header_column)
corrected_pvalues_table = calculate_corrected_pvalues(counts_table, pvalue_calculator, pvalue_corrector)
TableColumnCombiner.new([counts_table, corrected_pvalues_table]).combine.output($stdout)
