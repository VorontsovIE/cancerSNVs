require_relative 'genome_reader'

GENOME_FOLDER = File.absolute_path('../source_data/genome', __dir__)
EXONS_FILENAME = File.absolute_path('../source_data/exons.txt', __dir__)
MOTIF_NAMES_FILE = File.absolute_path('../source_data/motif_names.txt', __dir__)
GENE_INFOS = File.absolute_path('../source_data/hocomoco_genes_infos.csv', __dir__)
MOTIF_COLLECTION_FOLDER = File.absolute_path('../source_data/motif_collection/', __dir__)

KATAEGIS_COORDINATES_FILENAME = File.absolute_path('../source_data/AlexandrovEtAl/coordinates_of_kataegis.csv', __dir__)
SOMATIC_MUTATIONS_FOLDER = File.absolute_path('../source_data/AlexandrovEtAl/somatic_mutation_data', __dir__)
SAMPLE_INFOS_FILENAME = File.absolute_path('../source_data/AlexandrovEtAl/samples_summary.txt', __dir__)

ONE_BASED_INCLUSIVE = GenomeReader::CoordinateSystem::ONE_BASED_INCLUSIVE
ZERO_BASED_EXCLUSIVE = GenomeReader::CoordinateSystem::ZERO_BASED_EXCLUSIVE

# GENOME_READER = GenomeReader::DiskReader.new(
#   GENOME_FOLDER,
#   chromosome_file_by_name: ->(chr){ "chr#{chr}.plain" },
#   chromosome_name_matcher: /^chr(?<chromosome>\w+)\.plain$/
# )

GENOME_READER = GenomeReader::DiskReader.new(
  GENOME_FOLDER,
  chromosome_file_by_name: ->(chr){ "Homo_sapiens.GRCh37.75.dna_sm.chromosome.#{chr}.plain" },
  chromosome_name_matcher: /^Homo_sapiens\.GRCh37.75\.dna_sm\.chromosome\.(?<chromosome>\w+)\.plain$/
)
