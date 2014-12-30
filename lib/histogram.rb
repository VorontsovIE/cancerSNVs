# Use Rationals instead of floats!
# [from; to)
class Histogram
  attr_reader :from, :to, :step, :value_converter
  attr_reader :elements_total_in_range, :less_than, :greater_than

  # +original_from+ and +original_to+ are values from initial range
  # value_converter is used to convert both initial range boundaries and each value, which is put into a distribution
  def initialize(original_from,
                original_to,
                step,
                histogram_counts: nil,
                less_than: 0,
                greater_than: 0,
                elements_total_in_range: 0,
                &value_converter)
    @value_converter = block_given? ? value_converter : ->(val) { val }
    @original_from = original_from
    @original_to = original_to
    @step = step
    @from, @to = [@value_converter.call(original_from), @value_converter.call(original_to)].sort
    
    num_bins = (@from...@to).step(@step).size
    if histogram_counts
      unless histogram_counts.size == (@from...@to).step(@step).size
        raise "Histogram bins are incompatible with (from,to,step)-triple:\n" +
              "histogram_counts has #{histogram_counts.size} bins, while there are\n" +
              "#{bins.size} bins between #{from} and #{to} with step #{step}"
      end
      @histogram_counts = histogram_counts
    else
      @histogram_counts = Array.new(num_bins, 0)
    end

    @less_than = less_than
    @greater_than = greater_than
    @elements_total_in_range = elements_total_in_range
  end

  def elements_total
    @less_than + @greater_than + @elements_total_in_range
  end

  def clear
    @less_than = 0
    @greater_than = 0
    @elements_total_in_range = 0
    @histogram_counts = Array.new(num_bins, 0)
    self
  end

  def range
    @range ||= (from...to)
  end

  def num_bins
    @num_bins ||= (@from...@to).step(@step).size
  end

  def bins
    range.step(step).map{|left|
      (left...[left + step, to].min)
    }
  end

  # multiply each bin fold times
  def multiply(fold)
    Histogram.new(@original_from, @original_to, @step,
                  histogram_counts: @histogram_counts.map{|el| el * fold },
                  less_than: @less_than * fold,
                  greater_than: @greater_than * fold,
                  elements_total_in_range: @elements_total_in_range * fold,
                  &@value_converter)
  end

  def dup
    Histogram.new(@original_from, @original_to, @step,
                  histogram_counts: @histogram_counts.dup,
                  less_than: @less_than,
                  greater_than: @greater_than,
                  elements_total_in_range: @elements_total_in_range,
                  &@value_converter)
  end

  def each_bin(ignore_flanks: true)
    return enum_for(:each_bin, ignore_flanks: ignore_flanks)  unless block_given?

    yield (-Float::INFINITY...from), @less_than  unless ignore_flanks
    bins.each_with_index do |bin_range, index|
      yield bin_range, @histogram_counts[index]
    end
    yield (to..Float::INFINITY), @greater_than  unless ignore_flanks
  end

  def in_range?(value)
    range.include?(@value_converter.call(value))
  end

  def bin_index(value)
    value = @value_converter.call(value)
    return :less if value < from
    return :greater if value >= to
    ((value - from) / step).floor
  end
  private :bin_index

  def bin_count_for_value(value)
    index = bin_index(value)
    case index
    when :less
      @less_than
    when :greater
      @greater_than
    else
      @histogram_counts[index]
    end
  end

  def add_element(value)
    index = bin_index(value)
    # @elements_total += 1
    case index
    when :less
      @less_than += 1
      0 # element not added (more precisely, added out of range)
    when :greater
      @greater_than += 1
      0 # element not added (more precisely, added out of range)
    else
      @elements_total_in_range += 1
      @histogram_counts[index] += 1
      1 # element added
    end
  end
end
