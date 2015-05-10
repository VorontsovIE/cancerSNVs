Dir.glob(File.absolute_path('source_data', __dir__), '*.rake') do |fn|
  import fn
end

namespace 'source_data' do
  desc 'Download and prepare all source data. This includes somatic mutations, PerfectosAPE, genome assemblies and markup, motif collection'
  task prepare: ['Alexandrov:prepare', 'NikZainal:prepare', 'genome:prepare', 'exonic_markup:prepare', 'perfectosape']
end
