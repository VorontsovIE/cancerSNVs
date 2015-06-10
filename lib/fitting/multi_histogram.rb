require_relative 'histogram'
require_relative 'multi_histogram_fitter'

class MultiHistogram
  attr_reader :histograms, :block
  def initialize(histograms: {}, &block)
    raise 'Specify block which returns an empty (single, not indexed) histogram'  unless block_given?
    @histograms = histograms
    @block = block
  end

  def get_histogram(indices)
    if indices.size == 1
      @histograms[indices.first] ||= @block.call
    elsif indices.size > 1
      @histograms[indices.first] ||= MultiHistogram.new(&@block)
      @histograms[indices.first].get_histogram(indices.drop(1))
    else
      raise ArgumentError, 'indices array shouldn\'t be empty'
    end
  end

  def add_element(indices, value)
    get_histogram(indices).add_element(value)
  end

  def multiply(fold)
    multiplied_histograms = histograms.map{|k,v|
      [k, v.multiply(fold)]
    }.to_h
    MultiHistogram.new(histograms: multiplied_histograms, &@block)
  end

  def fitter(raise_on_missing: true)
    fitters = histograms.map{|k,v|
      [k, v.fitter]
    }.to_h
    MultiHistogramFitter.new(fitters, raise_on_missing: raise_on_missing)
  end

  def frequencies
    @histograms.sort.map{|index, inner_histogram| # inner_histogram can be multihistogram too
      [index, inner_histogram.frequencies]
    }.to_h
  end

  def to_s
    @histograms.sort.map{|index, inner_histogram| # inner_histogram can be multihistogram too
      [index, inner_histogram].map(&:to_s).join("\t")
    }.join("\n")
  end
end
