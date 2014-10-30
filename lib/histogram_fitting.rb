require_relative 'histogram'

class HistogramFitting
  attr_reader :goal_distribution, :current_distribution
  def initialize(sample_distribution)
    @goal_distribution = sample_distribution.empty_histogram
    @current_distribution = sample_distribution.empty_histogram
  end

  def add_element(object)
    goal_distribution.add_element(object)
  end

  def fit_element(object)
    count_current = current_distribution.bin_count_for_value(object)
    count_goal = goal_distribution.bin_count_for_value(object)
    if count_current < count_goal
      result = current_distribution.add_element(object)
      yield
      result # object added, either in range (then 0) or out of range (then 1)
    else
      0 # object not added because distribution already has necessary number of objects in an appropriate bin
    end
  end

  def goal_total
     @goal_distribution.elements_total_in_range
  end

  def current_total
     @current_distribution.elements_total_in_range
  end

  def print_discrepancies(msg = nil, output_stream = $stderr)
    unless current_distribution.elements_total == goal_distribution.elements_total
      output_stream.puts(msg)  if msg
      output_stream.puts "#{current_distribution.elements_total} < #{goal_distribution.elements_total}"
      (current_distribution.each_bin).zip(goal_distribution.each_bin).each do |(bin_range_1, bin_1_count), (bin_range_2, bin_2_count)|
        if bin_1_count != bin_2_count
          output_stream.puts(bin_range_1.round_to_s(3) + ": " + bin_1_count.to_s + " < " + bin_2_count.to_s)
        end
      end
    end
  end
end
