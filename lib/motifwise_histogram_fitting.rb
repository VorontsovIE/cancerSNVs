require_relative 'histogram_fitting'

class MotifHistogramFitter
  # histogram_smaple is a prototype of a histogram
  def initialize(motif_names, sample_distribution)
    @motif_names = motif_names
    @fitters = {}
    @motif_names.each do |motif_name|
      @fitters[motif_name] = HistogramFitting.new(sample_distribution)
    end
  end

  def add_element(motif_name, contexts, object)
    @fitters[motif_name].add_element(object)
  end
  def fit_element(motif_name, contexts, object, &block)
    @fitters[motif_name].fit_element(object, &block)
  end

  def goal_total
    @goal_total ||= @fitters.map{|motif_name, fitter| fitter.goal_total }.inject(0, &:+)
  end

  def current_total
    @fitters.map{|motif_name, fitter| fitter.current_total }.inject(0, &:+)
  end

  def print_discrepancies
    @motif_names.each do |motif_name|
      @fitters[motif_name][context_type].print_discrepancies("\n#{motif_name}", $stderr)
    end
  end
end
