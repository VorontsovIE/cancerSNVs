require 'rake'
require_relative 'lib/genome_reader'
require_relative 'lib/genome_markup'
require_relative 'lib/data_import/cancer_mutations_loading'

# Genome assembly, list of all exons and motif collection are huge
# thus it's possible to download them into any common place to be
# reusable across different projects.
# You can specify paths to common resources in corresponding files

# In project we use filesystem links to common resources so 
# one don't need to know system paths

# These can be safely redefined
module SystemPaths
  MotifCollection = File.absolute_path('./source_data/hocomoco', __dir__)
  MotifThresholds = File.absolute_path('./source_data/hocomoco-thresholds', __dir__)
  Genome          = File.absolute_path('./source_data/Ensembl-GRCh37.p13', __dir__)
  ExonicMarkup    = File.absolute_path('./source_data/hg19_exons(ensembl,GRCh37.p13).txt', __dir__)
  PerfectosAPE    = File.absolute_path('./source_data/ape.jar', __dir__)
end

# Attention! These should not be redefined!
module LocalPaths
  Root            = File.absolute_path(__dir__)

  MotifCollection = File.absolute_path('./source_data/motif_collection', __dir__)
  MotifThresholds = File.absolute_path('./source_data/motif_thresholds', __dir__)
  Genome          = File.absolute_path('./source_data/genome', __dir__)
  ExonicMarkup    = File.absolute_path('./source_data/exons.txt', __dir__)
  PerfectosAPE    = File.absolute_path('./ape.jar', __dir__)

  module Secondary
    AlexandrovData        = File.absolute_path('./source_data/AlexandrovEtAl', __dir__)
    CoordinatesOfKataegis = File.absolute_path('./source_data/AlexandrovEtAl/coordinates_of_kataegis.csv', __dir__)
    NikZainalSNVsOriginal = File.absolute_path('./source_data/SNV_infos_original.txt', __dir__)
    GeneInfos             = File.absolute_path('./source_data/hocomoco_genes_infos.csv', __dir__)
    MotifNames            = File.absolute_path('./source_data/motif_names.txt', __dir__) # To be removed

    SNVs                  = File.absolute_path('./results/SNVs', __dir__)
    Sequences             = File.absolute_path('./results/sequences', __dir__) # To be removed
    Chunks                = File.absolute_path('./results/sequence_chunks', __dir__)
    Sites                 = File.absolute_path('./results/sites')
    Fitting               = File.absolute_path('./results/fitted_sites')
    MotifStatistics       = File.absolute_path('./results/motif_statistics')

    module Alexandrov
      Mutations           = File.join(AlexandrovData, 'somatic_mutation_data')
      SamplesSummary      = File.join(AlexandrovData, 'samples_summary.txt')
    end

    NikZainalSNVs         = File.join(SNVs, 'SNV_infos_cancer.txt')
  end

  module Results
    module Alexandrov
      MutationTypesStatistics = File.join(Root, 'results/AlexandrovEtAl/mutation_types_statistics.txt')
      ContextStatistics       = File.join(Root, 'results/AlexandrovEtAl/mutation_contexts.txt')
    end
  end
end

module Configuration
  RandomGenomeFold = 100
  RandomShuffleFold = 100
  RandomGenomeSeeds = [13]
  RandomShuffleSeeds = [31]

  NumberOfCores = 8
  MemoryLimitOption = '' #'-Xmx512M'
  ExpandFlanksLength = 11
end

AlexandrovCancerTypes = FileList[File.join(LocalPaths::Secondary::Alexandrov::Mutations,'*')]
                          .select{|x| File.directory?(x) }
                          .pathmap('%f')
                          .to_a

AlexandrovWholeGenomeCancers = SampleInfo.each_in_file(LocalPaths::Secondary::Alexandrov::SamplesSummary)
                                .group_by(&:cancer_type)
                                .select{|cancer_type, samples| samples.any?(&:whole_genome?) }
                                .map{|cancer_type, samples| cancer_type }
                                .to_a.sort

ONE_BASED_INCLUSIVE = GenomeReader::CoordinateSystem::ONE_BASED_INCLUSIVE
ZERO_BASED_EXCLUSIVE = GenomeReader::CoordinateSystem::ZERO_BASED_EXCLUSIVE

# GENOME_READER = GenomeReader::DiskReader.new(
#   GENOME_FOLDER,
#   chromosome_file_by_name: ->(chr){ "chr#{chr}.plain" },
#   chromosome_name_matcher: /^chr(?<chromosome>\w+)\.plain$/
# )

GENOME_READER = GenomeReader::DiskReader.new(
  LocalPaths::Genome,
  chromosome_file_by_name: ->(chr){ "Homo_sapiens.GRCh37.75.dna_sm.chromosome.#{chr}.plain" },
  chromosome_name_matcher: /^Homo_sapiens\.GRCh37.75\.dna_sm\.chromosome\.(?<chromosome>\w+)\.plain$/
)

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
