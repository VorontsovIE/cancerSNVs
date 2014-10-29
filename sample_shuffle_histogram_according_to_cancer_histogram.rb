$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'histogram'

class Range
  def round_to_s(rate)
    from = self.begin.round(rate)
    to = self.end.round(rate)
    exclude_end? ? "#{from}...#{to}" : "#{from}..#{to}"
  end
end


def context_type(name_snp, context_types)
  context_types.select{|context_type, context_type_nameset| context_type_nameset.include?(name_snp) }
                .map{|context_type, context_type_nameset| context_type }
end

def load_context_types(mutations_filename, mutations_markup_filename)
  snps_splitted = File.readlines(mutations_filename).map{|el| el.chomp.split("\t")}
  cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
  tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
  not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
  any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }

  mut_types = File.readlines(mutations_markup_filename).drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
  intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
  promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
  regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type|
    intronic_mutation?(mut_type) || promoter_mutation?(mut_type)
  }

  { tpcpg: tpc_names & cpg_names & regulatory_mutation_names,
    cpg: cpg_names & regulatory_mutation_names,
    tpc: tpc_names & regulatory_mutation_names,
    not_cpg_tpc: not_cpg_tpc_names & regulatory_mutation_names#,
    # any_context: any_context_names & regulatory_mutation_names
  }
end

def each_regulatory_site(mutation_filename, regulatory_mutation_names, context_types)
  each_mutated_site_info(mutation_filename) do |mutated_site_info|
    motif_name = mutated_site_info.motif_name
    pvalue_1 = mutated_site_info.pvalue_1
    snp_name = mutated_site_info.normalized_snp_name
    next  unless pvalue_1 <= 0.0005
    next  unless regulatory_mutation_names.include?(snp_name)

    contexts = context_type(snp_name, context_types)
    yield mutated_site_info, motif_name, contexts, pvalue_1
  end
end

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

# range = from...to
create_histogram = ->{  Histogram.new(1e-7, 0.0005, 1.0){|pvalue| - Math.log2(pvalue) }  }

context_aware = true

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

context_types = load_context_types('source_data/SNPs.txt', 'source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt')

mutated_site_infos_cancer_filename = 'source_data/subsets/cancer_SNPs_regulatory_is_site.txt'
mutated_site_infos_shuffle_filename = 'source_data/subsets/shuffle_SNPs_regulatory_is_site.txt'
# mutated_site_infos_cancer_filename = 'source_data/cancer_SNPs.txt'
# mutated_site_infos_shuffle_filename = 'source_data/shuffle_SNPs.txt'

if context_aware
  fitters = ContextAwareMotifHistogramFitter.new(motif_names, context_types.keys, create_histogram.call)
else
  fitters = MotifHistogramFitter.new(motif_names, context_types.keys, create_histogram.call)
end


mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type|
  intronic_mutation?(mut_type) || promoter_mutation?(mut_type)
}


each_regulatory_site(mutated_site_infos_cancer_filename, regulatory_mutation_names, context_types) do |mutated_site_info, motif_name, contexts, pvalue_1|
  fitters.add_element(motif_name, contexts, pvalue_1)
end

$stderr.puts "Loaded #{fitters.goal_total}"

num_iteration = 0
$stderr.puts "Start shuffle reading"
loop do
  num_iteration += 1

  each_regulatory_site(mutated_site_infos_shuffle_filename, regulatory_mutation_names, context_types) do |mutated_site_info, motif_name, contexts, pvalue_1|
    fitters.fit_element(motif_name, contexts, pvalue_1) { puts mutated_site_info.line }
  end

  raise StopIteration  if fitters.current_total >= fitters.goal_total
  fitters.print_discrepancies

  $stderr.puts "Loaded #{fitters.current_total}"
  raise StopIteration # Now we run the only iteration
end
