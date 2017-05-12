require 'rake'
require_relative 'lib/genome_reader'
require_relative 'lib/genome_markup'
require_relative 'lib/data_import/cancer_mutations_loading'
require_relative 'lib/motif_family_recognizer'

# Genome assembly, list of all exons and motif collection are huge
# thus it's possible to download them into any common place to be
# reusable across different projects.
# You can specify paths to common resources in corresponding files

# In project we use filesystem links to common resources so
# one don't need to know system paths

# These can be safely redefined
module SystemPaths
  MotifCollection = File.absolute_path('./source_data/hocomoco', __dir__)
  MotifCollectionPCM = File.absolute_path('./source_data/hocomoco_pcm', __dir__)
  MotifThresholds = File.absolute_path('./source_data/hocomoco-thresholds', __dir__)
  Genome          = File.absolute_path('./source_data/Ensembl-GRCh37.p13', __dir__)
  ExonicMarkup    = File.absolute_path('./source_data/hg19_exons(ensembl,GRCh37.p13).txt', __dir__)
  PerfectosAPE    = File.absolute_path('./source_data/ape.jar', __dir__)
end

# Attention! These should not be redefined!
module LocalPaths
  Root            = File.absolute_path(__dir__)
  Results         = ENV['RESULTS_FOLDER'] || File.absolute_path('./results/', __dir__)

  MotifCollection = File.absolute_path('./source_data/motif_collection', __dir__)
  MotifCollectionPCM = File.absolute_path('./source_data/motif_collection_pcm', __dir__)
  MotifThresholds = File.absolute_path('./source_data/motif_thresholds', __dir__)
  Genome          = File.absolute_path('./source_data/genome', __dir__)
  ExonicMarkup    = File.absolute_path('./source_data/exons.txt', __dir__)
  PerfectosAPE    = File.absolute_path('./ape.jar', __dir__)

  module Secondary
    AlexandrovData        = File.absolute_path('./source_data/AlexandrovEtAl', __dir__)
    CoordinatesOfKataegis = File.absolute_path('./source_data/AlexandrovEtAl/coordinates_of_kataegis.csv', __dir__)
    MotifQualities        = File.absolute_path('./source_data/hocomoco_qualities.tsv', __dir__)
    HocomocoUniprots      = File.absolute_path('./source_data/HOCOMOCOv9_motifs2uniprot.txt', __dir__)
    UniprotDump           = File.absolute_path('./source_data/human_uniprot.txt', __dir__)
    MotifNames            = File.absolute_path('./source_data/motif_names.txt', __dir__) # To be removed

    SNVs                  = File.join(Results, 'SNVs')
    Sequences             = File.join(Results, 'sequences') # To be removed
    Chunks                = File.join(Results, 'sequence_chunks')
    Sites                 = File.join(Results, 'sites')
    Fitting               = File.join(Results, 'fitted_sites')
    MotifStatistics       = File.join(Results, 'motif_statistics')
    Slices                = File.join(Results, 'motif_statistics/slices')
    LogFolder             = File.join(Results, 'fitting_log')

    module Alexandrov
      Mutations           = File.join(AlexandrovData, 'somatic_mutation_data')
      SamplesSummary      = File.join(AlexandrovData, 'samples_summary.txt')

      MutationTypesStatisticsResults = File.join(Results, 'AlexandrovEtAl/mutation_types_statistics.txt')
      ContextStatisticsResults       = File.join(Results, 'AlexandrovEtAl/mutation_contexts.txt')
    end

  end

end

module Configuration
  # CorrectionMethod = 'BY' # Benjamini-Yekutieli
  CorrectionMethod = 'fdr' # FDR
  SignificanceCutoff = 0.05

  MotifQualities = 'A,B,C' # the only motif qualities to be taken into account

  FoldChangeCutoff = 4
  PvalueCutoff = 0.0005

  ExpandControlSetFold = 1 # one can artificially expand control set several time to quickly check whether expansion of control set will increase significance

  RandomGenomeSeeds = [] #[13,15,17]
  RandomShuffleSeeds = [] #[135,137,139]

  # Not using DHS profile data
  DHS_BED_FILES = { }

  # Using DHS profile data
  # DHS_BED_FILES = {:'Lung Adeno' => './source_data/dhs/LUNG_Roadmap_EncodeSet_merged.bed',
  #   :Breast => './source_data/dhs/BREAST_Roadmap_EncodeSet_merged.bed',
  #   :Liver => './source_data/dhs/LIVER_EncodeSet_merged.bed'}

  # Alexandrov works with the only seed
  AlexandrovRandomGenomeSeeds = 13
  AlexandrovRandomShuffleSeeds = 31

  module Alexandrov
    RandomGenomeDatasets = ['random_genome']
    RandomShuffleDatasets = ['random_shuffle']
    RandomDatasets = RandomGenomeDatasets + RandomShuffleDatasets
    Datasets = RandomDatasets + ['cancer']


    def self.contexts_by_cancer_type(cancer_type)
      [:any]
    end

    RandomGenomeFold = Hash.new(100).merge({:'Lung Adeno' => 10, :Breast => 30, :Liver => 20, :ALL => 1500, :AML => 3000, :'Pilocytic Astrocytoma' => 1000, :CLL => 200})
    RandomShuffleFold = Hash.new(100).merge({:'Lung Adeno' => 10, :Breast => 30, :Liver => 20, :ALL => 1500, :AML => 3000, :'Pilocytic Astrocytoma' => 1000, :CLL => 200})

    FittingFoldGenome = Hash.new(20).merge({:'Lung Adeno' => 2, :Breast => 6, :Liver => 4, :ALL => 300, :AML => 600, :'Pilocytic Astrocytoma' => 200, :CLL => 40})
    FittingFoldShuffle = Hash.new(20).merge({:'Lung Adeno' => 2, :Breast => 6, :Liver => 4, :ALL => 300, :AML => 600, :'Pilocytic Astrocytoma' => 200, :CLL => 40})

  # # DHS
  #  RandomGenomeFold = Hash.new(20).merge({:'Lung Adeno' => 10, :Breast => 5})
  #  RandomShuffleFold = Hash.new(20).merge({:'Lung Adeno' => 10, :Breast => 5})

  #  FittingFoldGenome = Hash.new(20).merge({:'Lung Adeno' => 2, :Breast => 1})
  #  FittingFoldShuffle = Hash.new(20).merge({:'Lung Adeno' => 2, :Breast => 1})
  end

  NumberOfCores = 24
  MemoryLimitOption = '-Xmx1G' # ''
  ExpandFlanksLength = 11

  def self.getAlexandrovWholeGenomeCancers
    @alexandrov_samples ||= begin
      if ENV['ALEXANDROV_SAMPLES']
        ENV['ALEXANDROV_SAMPLES'].split(',').map(&:to_sym).sort
      else
        SampleInfo.each_in_file(LocalPaths::Secondary::Alexandrov::SamplesSummary)
                  .group_by(&:cancer_type)
                  .select{|cancer_type, samples| samples.any?(&:whole_genome?) }
                  .map{|cancer_type, samples| cancer_type }
                  .to_a.sort
      end
    end
  end

  def self.sample_paths
    results = []
    results += Configuration.getAlexandrovWholeGenomeCancers.map{|cancer_type|
      ["#{cancer_type} (Alexandrov et al. sample)", File.join('Alexandrov', cancer_type.to_s)]
    }

    results.to_h
  end

  # part of pathname specifying sample with context for each cancer type for each experiment in each context
  def self.sample_with_context_paths
    results = []
    results += getAlexandrovWholeGenomeCancers.flat_map{|cancer_type|
      Alexandrov.contexts_by_cancer_type(cancer_type).map{|context|
        ["#{cancer_type} (Alexandrov et al. sample) in #{context} context", File.join('Alexandrov', cancer_type.to_s, context.to_s)]
      }
    }
    results.to_h
  end
end

ONE_BASED_INCLUSIVE = GenomeReader::CoordinateSystem::ONE_BASED_INCLUSIVE
ZERO_BASED_EXCLUSIVE = GenomeReader::CoordinateSystem::ZERO_BASED_EXCLUSIVE

GENOME_READER = GenomeReader::DiskReader.new(
  LocalPaths::Genome,
  chromosome_file_by_name: ->(chr){ "chr#{chr}.plain" },
  chromosome_name_matcher: /^chr(?<chromosome>\w+)\.plain$/
)

# GENOME_READER = GenomeReader::DiskReader.new(
#   LocalPaths::Genome,
#   chromosome_file_by_name: ->(chr){ "Homo_sapiens.GRCh37.75.dna_sm.chromosome.#{chr}.plain" },
#   chromosome_name_matcher: /^Homo_sapiens\.GRCh37.75\.dna_sm\.chromosome\.(?<chromosome>\w+)\.plain$/
# )

# One-based, inclusive
GENOME_MARKUP_LOADER = GenomeMarkupLoader.create(
  exonic_markup_filename: LocalPaths::ExonicMarkup,
  kataegis_coordinates_filename: LocalPaths::Secondary::CoordinatesOfKataegis,
  promoter_length_5_prime: 5000,
  promoter_length_3_prime: 500,
  kataegis_expansion_length: 1000
)

ALEXANDROV_MUTATIONS_LOADER = MutationsByCancerTypeLoader.create(
                                mutations_folder: LocalPaths::Secondary::Alexandrov::Mutations,
                                samples_summary_filename: LocalPaths::Secondary::Alexandrov::SamplesSummary
                              )

MOTIF_FAMILY_RECOGNIZERS = Hash.new{|hsh, deepness|
  hsh[deepness] = motif_family_recognizer_by_motif(
    deepness: deepness,
    tf_classification_filename: 'source_data/TFOntologies/TFClass_human.obo',
    uniprot_acs_by_id_filename: LocalPaths::Secondary::UniprotDump,
    uniprot_ids_by_motif_filename: LocalPaths::Secondary::HocomocoUniprots
  )
}
