# Use Rationals instead of floats!
# [from; to)
class Histogram
  attr_reader :from, :to, :step, :elements_total
  def range; @range ||= (from...to); end
  def bins; range.step(step).map{|left| (left...[left + step, to].min) }; end
  def each_bin
    return enum_for(:each_bin)  unless block_given?
    bins.each_with_index do |bin_range, index|
      # yield bin_range, @histogram[index]
      yield bin_range, @histogram_counts[index]
    end
  end
  def size; bins.size; end

  def initialize(from, to, step)
    @from, @to, @step = from, to, step
    # @histogram = Array.new(size) { [] }
    @histogram_counts = Array.new(size, 0)
    @elements_total = 0
  end

  def bin_index(value)
    return nil  unless range.include?(value)
    # raise ArgumentError, "value #{value} out of range [#{from}, #{to})"  unless range.include?(value)
    index = ((value - from) / step).floor
  end
  private :bin_index

  def bin_for_value(value)
    # @histogram[bin_index(value)]
    bin_index(value) && @histogram_counts[bin_index(value)]
  end

  def add_element(value)
    # bin_for_value(value) << value
    # bin_for_value(value) += 1
    return nil  unless bin_index(value)
    @histogram_counts[bin_index(value)] += 1
    @elements_total += 1
    value
  end
end
