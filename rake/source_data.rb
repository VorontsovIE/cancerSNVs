require_relative 'source_data/genome_assembly'
require_relative 'source_data/exonic_markup'
require_relative 'source_data/motif_collection'
require_relative 'source_data/perfectosape'

require_relative 'source_data/nik_zainal_data'
require_relative 'source_data/alexandrov_data'

namespace 'source_data' do
  desc 'Download and prepare all source data. This includes somatic mutations, PerfectosAPE, genome assemblies and markup, motif collection'
  task prepare: ['Alexandrov:prepare', 'NikZainal:prepare', 'genome:prepare', 'exonic_markup:prepare', 'perfectosape']
end
