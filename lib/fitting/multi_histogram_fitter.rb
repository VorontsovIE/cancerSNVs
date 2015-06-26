require_relative 'histogram_fitter'
require_relative 'multi_histogram'

class MultiHistogramFitter
  attr_reader :fitters
  # Constructor accepts a hash, where keys are motif names, values ar HistogramFitters
  def initialize(fitters, raise_on_missing: true)
    @fitters = fitters
    @raise_on_missing = raise_on_missing
  end

  def get_fitter(indices)
    raise ArgumentError, 'indices array shouldn\'t be empty'  if indices.empty?

    if indices.size == 1
      @fitters[indices.first] || nil
    elsif indices.size > 0
      @fitters[indices.first] ? @fitters[indices.first].get_fitter(indices.drop(1)) : nil
    end
  end

  def element_can_be_added?(indices, object)
    fitter = get_fitter(indices)
    fitter && fitter.element_can_be_added?(object)
  end

  def fit_element(indices, object, &block)
    fitter = get_fitter(indices)
    if fitter
      fitter.fit_element(object, &block)
    elsif @raise_on_missing # otherwise do nothing
      raise 'An attept to fit an element which never was in an original distribution'
    end
  end

  def goal_total
    @goal_total ||= @fitters.each_value.map(&:goal_total).inject(0, &:+)
  end

  def current_total
    @fitters.each_value.map(&:current_total).inject(0, &:+)
  end

  def num_underfitted
    goal_total - current_total
  end

  def underfitting_percentage
    100.0 * num_underfitted / goal_total
  end

  def fitted?
    current_total >= goal_total
  end


  def print_discrepancies(output_stream: $stderr, note: nil, level_mark: "")
    return  if fitted?

    output_stream.puts "#{level_mark}#{note}"
    output_stream.puts "#{level_mark}#{num_underfitted} underfitted (#{underfitting_percentage}%)"
    @fitters.sort.each do |index, fitter|
      fitter.print_discrepancies(output_stream: output_stream, note: index, level_mark: "#{level_mark}\t")
    end
  end
end
