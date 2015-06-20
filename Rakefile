require_relative 'experiment_configuration'

AlexandrovWholeGenomeCancers = Configuration.getAlexandrovWholeGenomeCancers
YeastApobecSamples = Configuration.getYeastApobecSamples

task :load_genome_markup do
  GENOME_MARKUP ||= GENOME_MARKUP_LOADER.load_markup
end

import *FileList['rakelib/processing/*.rake']
