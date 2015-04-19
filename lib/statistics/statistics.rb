# imported from EpigeneDB
class Statistics
  attr_reader :values

  def initialize(values)
    @values = values
  end

  def size
    values.size
  end

  def sorted_values
    values.sort
  end

  def median
    return nil  if values.empty?
    sorted_values[size / 2]
  end

  def mean
    return nil  if values.empty?
    values.inject(0.0, &:+) / size
  end

  def mean_square
    return nil  if values.empty?
    values.map{|x| x**2}.inject(0.0, &:+) / size
  end

  def variance
    return nil  if values.empty?
    result = mean_square - mean**2
    result >= 0 ? result : 0 # floating point bug
  end

  def stddev
    return nil  if size <= 1
    Math.sqrt((size / (size - 1.0)) * variance)
  end

  def min
    return nil  if values.empty?
    values.min
  end

  def max
    return nil  if values.empty?
    values.max
  end

  def relative_stddev
    return nil  unless mean && stddev && mean != 0
    stddev / mean
  end
end
