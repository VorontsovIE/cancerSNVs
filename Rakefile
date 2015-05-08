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
    NikZainalSNVs         = File.absolute_path('./source_data/SNV_infos_original.txt', __dir__)
    GeneInfos             = File.absolute_path('./source_data/hocomoco_genes_infos.csv', __dir__)
    MotifNames            = File.absolute_path('./source_data/motif_names.txt', __dir__) # To be removed

    SNVs                  = File.absolute_path('./results/SNVs', __dir__) # To be removed
    Sequences             = File.absolute_path('./results/sequences', __dir__)
    Chunks                = File.absolute_path('./results/sequence_chunks', __dir__)
    Sites                 = File.absolute_path('./results/sites')
    Fitting               = File.absolute_path('./results/fitted_sites')
    MotifStatistics       = File.absolute_path('./results/motif_statistics')
  end
end

require_relative 'rake/source_data'
require_relative 'rake/preparations'
