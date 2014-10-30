require_relative 'histogram_fitting'

class ContextAwareMotifHistogramFitter
  # histogram_smaple is a prototype of a histogram
  def initialize(motif_names, context_types, sample_distribution)
    @motif_names = motif_names
    @context_types = context_types
    @fitters = {}
    @motif_names.each do |motif_name|
      @fitters[motif_name] = {}
      @context_types.each do |context_type|
        @fitters[motif_name][context_type] = HistogramFitting.new(sample_distribution)
      end
    end
  end

  def add_element(motif_name, context_types, object)
    context_types.each do |context_type|
      @fitters[motif_name][context_type].add_element(object)
    end
  end
  def fit_element(motif_name, context_types, object, &block)
    context_types.each do |context_type|
      @fitters[motif_name][context_type].fit_element(object, &block)
    end
  end

  def goal_total
    @goal_total ||= @fitters.map{|motif_name, motif_fitters|
      motif_fitters.map{|context_type, fitter|
        fitter.goal_total
      }.inject(0, &:+)
    }.inject(0, &:+)
  end

  def current_total
    @fitters.map{|motif_name, motif_fitters|
      motif_fitters.map{|context_type, fitter|
        fitter.current_total
      }.inject(0, &:+)
    }.inject(0, &:+)
  end

  def print_discrepancies
    @motif_names.each do |motif_name|
      @context_types.each do |context_type|
        @fitters[motif_name][context_type].print_discrepancies("\n#{motif_name},#{context_type}", $stderr)
      end
    end
  end
end
