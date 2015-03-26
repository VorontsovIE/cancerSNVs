require 'open3'
require 'rubystats'
require_relative 'table'

class PvalueCalculator
  attr_reader :class_counts
  def initialize(class_counts: :two_classes)
    raise 'Class counts correction can be either :two_classes or :class_and_total'  unless [:two_classes, :class_and_total].include?(class_counts)
    @class_counts = class_counts
  end

  # either  a_1, b_1, a_2, b_2
  # or      a_1, b_1, total_a, total_b
  # depending on class_counts
  # unclassified are not included in total, they are additional
  def calculate(values, unclassified_a: 0, unclassified_b: 0)
    case class_counts
    when :two_classes
      a_1, b_1, a_2, b_2 = *values
    when :class_and_total
      a_1, b_1, total_a, total_b = *values
      a_2 = total_a - a_1
      b_2 = total_b - b_1
    end

    fet = Rubystats::FishersExactTest.new
    # Sorry, slow method. I don't want to make assumptions about significance
    (0..unclassified_a).flat_map{|a_1_add|
      a_2_add = unclassified_a - a_1_add

      (0..unclassified_b).map{|b_1_add|
        b_2_add = unclassified_b - b_1_add

        fet.calculate(a_1 + a_1_add, a_2 + a_2_add,
                      b_1 + b_1_add, b_2 + b_2_add)[:twotail]
      }
    }.max
  end
end

class PvalueCorrector
  attr_reader :correction_method
  def initialize(correction_method)
    @correction_method = correction_method
  end

  def correct(pvalues)
    script_path = File.absolute_path('correct_significance.r', __dir__) # runned with Rscript (shebang syntax involved)
    Open3.popen2("#{script_path} #{correction_method}"){|fread, fwrite|
      fread.puts pvalues
      fread.close
      fwrite.readlines.map(&:strip).map(&:to_f)
    }
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