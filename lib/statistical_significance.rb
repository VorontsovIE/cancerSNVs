require 'tempfile'
require 'rubystats'

def calculate_fisher_test(disrupted_shuffle, disrupted_cancer, total_shuffle, total_cancer)
  Rubystats::FishersExactTest.new.calculate(
                                            disrupted_shuffle,
                                            total_shuffle - disrupted_shuffle,
                                            disrupted_cancer,
                                            total_cancer - disrupted_cancer
                                          )[:twotail]
end


def with_temp_file(filename, &block)
  temp_file = Tempfile.new(filename)
  yield temp_file
ensure
  temp_file.close
  temp_file.unlink
end


def calculate_holms_correction(pvalues)
  with_temp_file('uncorrected_pvalues.txt') do |uncorrected_pvalues_file|
    pvalues.each do |pvalue| 
      uncorrected_pvalues_file.puts(pvalue)
    end
    uncorrected_pvalues_file.close

    with_temp_file('corrected_pvalues.txt') do |corrected_pvalues_file|
      corrected_pvalues_file.close

      with_temp_file('correction_script.r') do |correction_script_file|
        correction_script_file.puts('dataT <- read.table("' + uncorrected_pvalues_file.path + '")')
        correction_script_file.puts('values <- data.frame(dataT)')
        correction_script_file.puts('newValues <- p.adjust(values$V1, "holm")')
        correction_script_file.puts('write.table(newValues, "'+ corrected_pvalues_file.path + '", FALSE, FALSE, "\t", "\n", "NA", ".", FALSE, FALSE)')
        correction_script_file.close

        system("/usr/bin/env R --slave -f #{correction_script_file.path}")  
      end

      corrected_pvalues_file.open
      corrected_pvalues_file.readlines.map(&:strip).map(&:to_f)
    end
  end
end
