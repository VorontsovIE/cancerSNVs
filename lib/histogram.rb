# Use Rationals instead of floats!
# [from; to)
class Histogram
  attr_reader :from, :to, :step, :value_converter, :elements_total, :elements_total_in_range, :less_than, :greater_than

  # from and to are values from initial range
  # value_converter is used to convert both initial range boundaries and each value, which is put into a distribution
  def initialize(from, to, step, &value_converter)
    if block_given?
      @value_converter = value_converter
    else
      @value_converter = ->(val) { val }
    end
    @original_from = from
    @original_to = to
    @from, @to = [@value_converter.call(from), @value_converter.call(to)].sort
    @step = step
    @histogram_counts = Array.new(bins.size, 0)

    @less_than = 0
    @greater_than = 0

    @elements_total = 0
    @elements_total_in_range = 0
  end

  def empty_histogram
    Histogram.new(@original_from, @original_to, @step, &@value_converter)
  end

  def range
    @range ||= (from...to)
  end

  def bins
    range.step(step).map{|left|
      (left...[left + step, to].min)
    }
  end

  def each_bin(ignore_flanks: false)
    return enum_for(:each_bin, ignore_flanks: ignore_flanks)  unless block_given?

    yield (-Float::INFINITY...from), @less_than  unless ignore_flanks
    bins.each_with_index do |bin_range, index|
      yield bin_range, @histogram_counts[index]
    end
    yield (to..Float::INFINITY), @greater_than  unless ignore_flanks
  end

  # def size; bins.size; end


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
    @elements_total += 1
    case index
    when :less
      @less_than += 1
      0 # element not added (more precisely, added out of range)
    when :greater
      @greater_than += 1
      0 # element not added (more precisely, added out of range)
    else
      @elements_total_in_range +=1
      @histogram_counts[index] += 1
      1 # element added
    end
  end
end
