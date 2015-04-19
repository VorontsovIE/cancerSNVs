require 'forwardable'
require_relative 'statistics/fisher_table'
require_relative 'statistics/statistical_significance'

def to_float(str)
  str && str.to_f
end

class MotifStatistics
  attr_reader :motif
  attr_reader :disruption_table, :emergence_table, :core_flank_table
  attr_reader :random_unclassified
  attr_reader :gene, :quality, :official_gene_name

  def initialize( motif:,
                  disruption_table:, emergence_table:, core_flank_table:,
                  random_unclassified:,
                  gene:, quality:, official_gene_name:
                )
    @motif = motif
    @disruption_table = disruption_table
    @emergence_table = emergence_table
    @core_flank_table = core_flank_table
    @random_unclassified = random_unclassified

    @gene = gene
    @quality = quality
    @official_gene_name = official_gene_name
  end

  #####
  def underfitted; random_unclassified; end
  def cancer_total_before_substitution; disruption_table.class_a_total ; end
  def random_total_before_substitution; disruption_table.class_b_total ; end

  def cancer_disrupted; disruption_table.class_a_positive ; end
  def random_disrupted; disruption_table.class_b_positive ; end

  def cancer_core; core_flank_table.class_a_positive; end
  def cancer_flank; core_flank_table.class_a_negative; end
  def random_core; core_flank_table.class_b_positive; end
  def random_flank; core_flank_table.class_b_negative; end

  #####
  def cancer_total_after_substitution; emergence_table.class_a_total ; end
  def random_total_after_substitution; emergence_table.class_b_total ; end

  def cancer_emerged; emergence_table.class_a_positive ; end
  def random_emerged; emergence_table.class_b_positive ; end

  #####
  def cancer_disruption_rate; disruption_table.class_a_positive_rate; end
  def random_disruption_rate; disruption_table.class_b_positive_rate; end
  def cancer_to_random_disruption_ratio; disruption_table.a_to_b_positive_rate_ratio; end

  def cancer_emergence_rate; emergence_table.class_a_positive_rate; end
  def random_emergence_rate; emergence_table.class_b_positive_rate; end
  def cancer_to_random_emergence_ratio; emergence_table.a_to_b_positive_rate_ratio; end

  #####
  def cancer_core_rate; core_flank_table.class_a_positive_rate; end
  def random_core_rate; core_flank_table.class_b_positive_rate; end
  def cancer_to_random_core_ratio; core_flank_table.a_to_b_positive_rate_ratio; end

  #####
  def disruption_significance; disruption_table.significance; end
  def disruption_significance_fitting_aware; disruption_table.significance(unclassified_b: random_unclassified); end

  def emergence_significance; emergence_table.significance; end
  def emergence_significance_fitting_aware; emergence_table.significance(unclassified_b: random_unclassified); end

  def core_flank_significance; core_flank_table.significance; end
  def core_flank_significance_fitting_aware; core_flank_table.significance(unclassified_b: random_unclassified); end
end


class MotifCollectionStatistics
  attr_reader :collection
  def initialize(collection, pvalue_corrector: IdlePvalueCorrector.new)
    unless collection.all?{|motif_statistics| motif_statistics.is_a?(MotifStatistics) }
      raise "MotifCollectionStatistics wraps MotifStatistics"
    end
    @collection = collection.map{|motif_statistics| MotifStatisticsCorrected.new(motif_statistics) }

    disruption_significances = pvalue_corrector.correct( @collection.map(&:disruption_significance_uncorrected) )
    emergence_significances = pvalue_corrector.correct( @collection.map(&:emergence_significance_uncorrected) )
    core_flank_significances = pvalue_corrector.correct( @collection.map(&:core_flank_significance_uncorrected) )

    disruption_significances_fitting_aware = pvalue_corrector.correct( @collection.map(&:disruption_significance_uncorrected_fitting_aware) )
    emergence_significances_fitting_aware = pvalue_corrector.correct( @collection.map(&:emergence_significance_uncorrected_fitting_aware) )
    core_flank_significances_fitting_aware = pvalue_corrector.correct( @collection.map(&:core_flank_significance_uncorrected_fitting_aware) )

    @collection.each_with_index do |motif_statistics, ind|
      motif_statistics.disruption_significance = disruption_significances[ind]
      motif_statistics.emergence_significance = emergence_significances[ind]
      motif_statistics.disruption_significance_fitting_aware = disruption_significances_fitting_aware[ind]
      motif_statistics.emergence_significance_fitting_aware = emergence_significances_fitting_aware[ind]
      motif_statistics.core_flank_significance = core_flank_significances[ind]
      motif_statistics.core_flank_significance_fitting_aware = core_flank_significances_fitting_aware[ind]
    end
  end

  def self.table_header
    COLUMN_ORDER.map{|column_id| COLUMN_NAMES[column_id] }.join("\t")
  end

  def to_s
    ([MotifCollectionStatistics.table_header] + collection.map(&:to_s)).join("\n")
  end


  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.readline # skip header
      each_in_stream(f, &block)
    end
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.map{|line| MotifStatisticsCorrected.from_string(line) }.each(&block)
  end

  class MotifStatisticsCorrected
    extend Forwardable
    def_delegators :@motif_statistics_uncorrected, :motif, :gene, :quality, :official_gene_name, :underfitted,
                                                   :cancer_total_before_substitution, :random_total_before_substitution,
                                                   :cancer_disrupted, :random_disrupted,
                                                   :cancer_total_after_substitution, :random_total_after_substitution,
                                                   :cancer_emerged, :random_emerged,
                                                   :cancer_core, :cancer_flank, :random_core, :random_flank,
                                                   :cancer_disruption_rate, :random_disruption_rate, :cancer_to_random_disruption_ratio,
                                                   :cancer_emergence_rate, :random_emergence_rate, :cancer_to_random_emergence_ratio,
                                                   :cancer_core_rate, :random_core_rate, :cancer_to_random_core_ratio

    attr_reader :motif_statistics_uncorrected
    def initialize(motif_statistics_uncorrected, &block)
      @motif_statistics_uncorrected = motif_statistics_uncorrected
      instance_eval(&block)  if block_given?
    end

    attr_accessor :disruption_significance, :disruption_significance_fitting_aware,
                  :emergence_significance, :emergence_significance_fitting_aware,
                  :core_flank_significance, :core_flank_significance_fitting_aware
    attr_writer :disruption_significance_uncorrected, :disruption_significance_uncorrected_fitting_aware,
                :emergence_significance_uncorrected, :emergence_significance_uncorrected_fitting_aware,
                :core_flank_significance_uncorrected, :core_flank_significance_uncorrected_fitting_aware
    def disruption_significance_uncorrected
      @disruption_significance_uncorrected ||= motif_statistics_uncorrected.disruption_significance
    end

    def disruption_significance_uncorrected_fitting_aware
      @disruption_significance_uncorrected_fitting_aware ||= motif_statistics_uncorrected.disruption_significance_fitting_aware
    end

    def emergence_significance_uncorrected
      @emergence_significance_uncorrected ||= motif_statistics_uncorrected.emergence_significance
    end

    def emergence_significance_uncorrected_fitting_aware
      @emergence_significance_uncorrected_fitting_aware ||= motif_statistics_uncorrected.emergence_significance_fitting_aware
    end

    def core_flank_significance_uncorrected
      @core_flank_significance_uncorrected ||= motif_statistics_uncorrected.core_flank_significance
    end

    def core_flank_significance_uncorrected_fitting_aware
      @core_flank_significance_uncorrected_fitting_aware ||= motif_statistics_uncorrected.core_flank_significance_fitting_aware
    end

    def to_s
      COLUMN_ORDER.map{|meth| send(meth) }.join("\t")
    end

    def self.from_string(line)
      motif, official_gene_name, quality,
          _cancer_to_random_disruption_ratio, disruption_significance_fitting_aware,
          _cancer_to_random_emergence_ratio, emergence_significance_fitting_aware,
          _cancer_to_random_core_ratio, core_flank_significance_fitting_aware,
          random_disrupted, random_emerged, random_total_before_substitution, random_total_after_substitution, random_core, random_flank,
          cancer_disrupted, cancer_emerged, cancer_total_before_substitution, cancer_total_after_substitution, cancer_core, cancer_flank,
          _cancer_disruption_rate, _cancer_emergence_rate, _cancer_core_rate, _random_disruption_rate, _random_emergence_rate, _random_core_rate,
          disruption_significance_uncorrected, emergence_significance_uncorrected, core_flank_significance_uncorrected,
          gene, underfitted,
          disruption_significance, emergence_significance, core_flank_significance,
          disruption_significance_uncorrected_fitting_aware, emergence_significance_uncorrected_fitting_aware, core_flank_significance_uncorrected_fitting_aware = line.chomp.split("\t")

      motif_statistics = MotifStatistics.new(
        motif: motif.to_sym,
        disruption_table: FisherTable.by_class_and_total(
          class_a_total: cancer_total_before_substitution.to_i, class_a_positive: cancer_disrupted.to_i,
          class_b_total: random_total_before_substitution.to_i,  class_b_positive: random_disrupted.to_i,
        ),
        emergence_table: FisherTable.by_class_and_total(
          class_a_total: cancer_total_after_substitution.to_i,  class_a_positive: cancer_emerged.to_i,
          class_b_total: random_total_after_substitution.to_i,  class_b_positive: random_emerged.to_i,
        ),
        core_flank_table: FisherTable.by_two_classes(
          class_a_positive: cancer_core.to_i,  class_a_negative: cancer_flank.to_i,
          class_b_positive: random_core.to_i,  class_b_negative: random_flank.to_i,
        ),
        random_unclassified: underfitted.to_i,
        gene: gene, quality: quality.upcase.to_sym, official_gene_name: official_gene_name
      )

      result = MotifStatisticsCorrected.new(motif_statistics) do |s|
        s.disruption_significance = to_float(disruption_significance)
        s.disruption_significance_fitting_aware = to_float(disruption_significance_fitting_aware)
        s.emergence_significance = to_float(emergence_significance)
        s.emergence_significance_fitting_aware = to_float(emergence_significance_fitting_aware)
        s.core_flank_significance = to_float(core_flank_significance)
        s.core_flank_significance_fitting_aware = to_float(core_flank_significance_fitting_aware)

        s.disruption_significance_uncorrected = to_float(disruption_significance_uncorrected) # it can be recalculated, but it takes time
        s.emergence_significance_uncorrected = to_float(emergence_significance_uncorrected)
        s.disruption_significance_uncorrected_fitting_aware = to_float(disruption_significance_uncorrected_fitting_aware)
        s.emergence_significance_uncorrected_fitting_aware = to_float(emergence_significance_uncorrected_fitting_aware)
        s.core_flank_significance_uncorrected = to_float(core_flank_significance_uncorrected)
        s.core_flank_significance_uncorrected_fitting_aware = to_float(core_flank_significance_uncorrected_fitting_aware)
      end
    end
  end

  COLUMN_NAMES = {
    motif: 'Motif',
    quality: 'Motif quality',
    gene: 'TF gene',
    official_gene_name: 'TF gene (official name)',
    random_disrupted: 'Random disrupted',
    random_emerged: 'Random emerged',
    random_total_before_substitution: 'Random total sites before substitution',
    random_total_after_substitution: 'Random total sites after substitution',
    cancer_disrupted: 'Cancer disrupted',
    cancer_emerged: 'Cancer emerged',
    cancer_total_before_substitution: 'Cancer total sites before substitution',
    cancer_total_after_substitution: 'Cancer total sites after substitution',

    cancer_core: 'Cancer substitutions in motif core',
    cancer_flank: 'Cancer substitutions in motif flank',
    random_core: 'Random substitutions in motif core',
    random_flank: 'Random substitutions in motif flank',
    cancer_core_rate: 'Fraction of cancer substitutions hitted core of motif',
    random_core_rate: 'Fraction of random substitutions hitted core of motif',
    cancer_to_random_core_ratio: 'Cancer to random substitutions in core of motif ratio',
    core_flank_significance_uncorrected: 'Significance of difference in core-flank rate (not corrected on multiple comparison, no correction for underfitted sites)',
    core_flank_significance_uncorrected_fitting_aware: 'Significance of difference in core-flank rate (no correction for multiple comparison, fitting aware)',
    core_flank_significance: 'Significance of difference in core-flank rate (corrected on multiple comparison, no correction for underfitted sites)',
    core_flank_significance_fitting_aware: 'Significance of difference in core-flank rate (corrected on multiple comparison, fitting aware)',

    cancer_disruption_rate: 'Cancer disruption rate',
    cancer_emergence_rate: 'Cancer emergence rate',
    random_disruption_rate: 'Random disruption rate',
    random_emergence_rate: 'Random emergence rate',
    cancer_to_random_disruption_ratio: 'Cancer to random disruption rate ratio',
    cancer_to_random_emergence_ratio: 'Cancer to random emergence rate ratio',
    disruption_significance_uncorrected: 'Significance of difference in disruption rate (not corrected)',
    emergence_significance_uncorrected: 'Significance of difference in emergence rate (not corrected)',
    disruption_significance: 'Significance of difference in disruption rate (corrected on multiple comparison)',
    emergence_significance: 'Significance of difference in emergence rate (corrected on multiple comparison)',
    underfitted: 'Number of sites underfitted',
    disruption_significance_uncorrected_fitting_aware: 'Significance of difference in disruption rate (not corrected, fitting aware)',
    emergence_significance_uncorrected_fitting_aware: 'Significance of difference in emergence rate (not corrected, fitting aware)',
    disruption_significance_fitting_aware: 'Significance of difference in disruption rate (corrected on multiple comparison, fitting aware)',
    emergence_significance_fitting_aware: 'Significance of difference in emergence rate (corrected on multiple comparison, fitting aware)',
  }

  COLUMN_ORDER = [:motif, :official_gene_name, :quality,
                  :cancer_to_random_disruption_ratio, :disruption_significance_fitting_aware,
                  :cancer_to_random_emergence_ratio, :emergence_significance_fitting_aware,
                  :cancer_to_random_core_ratio, :core_flank_significance_fitting_aware,
                  :random_disrupted, :random_emerged, :random_total_before_substitution, :random_total_after_substitution, :random_core, :random_flank,
                  :cancer_disrupted, :cancer_emerged, :cancer_total_before_substitution, :cancer_total_after_substitution, :cancer_core, :cancer_flank,
                  :cancer_disruption_rate, :cancer_emergence_rate, :cancer_core_rate, :random_disruption_rate, :random_emergence_rate, :random_core_rate,
                  :disruption_significance_uncorrected, :emergence_significance_uncorrected, :core_flank_significance_uncorrected,
                  :gene, :underfitted,
                  :disruption_significance, :emergence_significance, :core_flank_significance,
                  :disruption_significance_uncorrected_fitting_aware, :emergence_significance_uncorrected_fitting_aware, :core_flank_significance_uncorrected_fitting_aware,
                ]
end
