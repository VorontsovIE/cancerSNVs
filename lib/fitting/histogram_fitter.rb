require_relative 'histogram'

class HistogramFitter
  attr_reader :goal_distribution, :current_distribution
  def initialize(goal_distribution)
    @goal_distribution = goal_distribution
    @current_distribution = goal_distribution.dup
    @current_distribution.clear
  end

  def element_can_be_added?(object)
    count_current = current_distribution.bin_count_for_value(object)
    count_goal = goal_distribution.bin_count_for_value(object)
    count_current < count_goal
  end

  def fit_element(object)
    if element_can_be_added?(object)
      current_distribution.add_element(object)
      yield
    end
  end

  def goal_total
     @goal_distribution.elements_total_in_range
  end

  def current_total
     @current_distribution.elements_total_in_range
  end

  def fitted?
    current_total >= goal_total
  end

  def num_underfitted
    goal_total - current_total
  end

  def underfitting_percentage
    100.0 * num_underfitted / goal_total
  end

  def each_bin(ignore_flanks: true)
    return enum_for(:each_bin, ignore_flanks: ignore_flanks)  unless block_given?
    iter_current = current_distribution.each_bin(ignore_flanks: ignore_flanks)
    iter_goal = goal_distribution.each_bin(ignore_flanks: ignore_flanks)
    iter_current.zip(iter_goal).each do |(bin_range, current_bin_count), (_same_bin_range, goal_bin_count)|
      yield bin_range, current_bin_count, goal_bin_count
    end
  end

  def print_discrepancies(output_stream: $stderr, note: nil, level_mark: "")
    return  if fitted?

    output_stream.puts "#{level_mark}#{note}"
    output_stream.puts "#{level_mark}#{num_underfitted} underfitted (#{underfitting_percentage}%)"
    each_bin do |bin_range, current_bin_count, goal_bin_count|
      if current_bin_count != goal_bin_count
        range_repr = range_formatting(bin_range, rate: 3)
        output_stream.puts "#{level_mark}#{range_repr}: found #{current_bin_count} from #{goal_bin_count}"
      end
    end
    # output_stream.puts "#{level_mark}"
  end
end
