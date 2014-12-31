require 'tempfile'
require 'rubystats'
require_relative 'table'

def with_temp_file(filename, &block)
  temp_file = Tempfile.new(filename)
  yield temp_file
ensure
  temp_file.close
  temp_file.unlink
end

class PvalueCalculator
  attr_reader :class_counts
  def initialize(class_counts: :two_classes)
    raise 'Class counts correction can be either :two_classes or :class_and_total'  unless [:two_classes, :class_and_total].include?(class_counts)
    @class_counts = class_counts
  end

  # either  n_11,n_21,n_12, n_22
  # or      n_11, n_21, total_1, total_2
  # depending on class_counts
  def calculate(values)
    case class_counts
    when :two_classes
      n_11, n_21, n_12, n_22 = *values
      Rubystats::FishersExactTest.new.calculate(n_11, n_12, n_21, n_22)[:twotail]
    when :class_and_total
      n_11, n_21, total_1, total_2 = *values
      Rubystats::FishersExactTest.new.calculate(n_11, total_1 - n_11, n_21, total_2 - n_21)[:twotail]
    end
  end
end

class HolmsPvalueCorrector
  def correct(pvalues)
    with_temp_file('uncorrected_pvalues.txt') do |uncorrected_pvalues_file|
      pvalues.each do |pvalue|
        uncorrected_pvalues_file.puts(pvalue)
      end
      uncorrected_pvalues_file.close

      with_temp_file('corrected_pvalues.txt') do |corrected_pvalues_file|
        corrected_pvalues_file.close

        run_correction_script(uncorrected_pvalues_file.path, corrected_pvalues_file.path)

        corrected_pvalues_file.open
        corrected_pvalues_file.readlines.map(&:strip).map(&:to_f)
      end
    end
  end

  def correct_hash(pvalues_hash)
    names = []
    pvalues = []
    pvalues_hash.each{|name, pvalue|
      names << name
      pvalues << pvalue
    }
    names.zip(correct(pvalues)).to_h
  end

  def run_correction_script(from_file, to_file)
    with_temp_file('correction_script.r') do |correction_script_file|
      correction_script_file.puts('dataT <- read.table("' + from_file + '")')
      correction_script_file.puts('values <- data.frame(dataT)')
      correction_script_file.puts('newValues <- p.adjust(values$V1, "holm")')
      correction_script_file.puts('write.table(newValues, "'+ to_file + '", FALSE, FALSE, "\t", "\n", "NA", ".", FALSE, FALSE)')
      correction_script_file.close

      system("/usr/bin/env R --slave -f #{correction_script_file.path}")
    end
  end
  private :run_correction_script
end

class IdlePvalueCorrector
  def correct(pvalues)
    pvalues
  end
end

def calculate_corrected_pvalues(counts_table, pvalue_calculator, pvalue_corrector, header: [:pvalue, 'P-value'])
  pvalues = counts_table.each.map do |line|
    class_counts = line.map(&:to_i)
    raise "Bad input, should be 4 class count columns. See line:\n#{line.inspect}"  unless class_counts.size == 4
    pvalue_calculator.calculate(class_counts)
  end
  corrected_pvalues = pvalue_corrector.correct(pvalues)
  Table.new(corrected_pvalues.map{|value| [value] }, header_row: [header], header_column: nil)
end