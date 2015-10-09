require_relative '../../lib/calculate_contexts'

desc 'Collect sample statistics (all SNVs, not only regulatory)'
task :sample_statistics_all_SNVs do
  header = ['Sample', 'Genomic control fold', 'Shuffle control fold', 'Genomic control fitting fold', 'Shuffle control fitting fold', 'Mutations total']
  matrix = [ [*header, *possible_contexts] ]
  short_matrix = [ [*header, *possible_short_contexts] ]
  Configuration.getAlexandrovWholeGenomeCancers.each{|cancer_type|
    context_counts = context_counts_in_file(File.join(LocalPaths::Results, 'AllSNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt'))
    short_context_counts = short_context_counts_in_file(File.join(LocalPaths::Results, 'AllSNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt'))
    cancer_infos = ["#{cancer_type} (Alexandrov et al. sample)",
                    Configuration::Alexandrov::RandomGenomeFold[cancer_type],
                    Configuration::Alexandrov::RandomShuffleFold[cancer_type],
                    Configuration::Alexandrov::FittingFoldGenome[cancer_type],
                    Configuration::Alexandrov::FittingFoldShuffle[cancer_type]]

    matrix << [*cancer_infos,
              context_counts.each_value.inject(0, &:+),
              *possible_contexts.map{|context| context_counts[context] }]
    short_matrix << [*cancer_infos,
                    short_context_counts.each_value.inject(0, &:+),
                    *possible_short_contexts.map{|context| short_context_counts[context] }]
  }

  Configuration.getYeastApobecSamples.each{|cancer_type|
    context_counts = context_counts_in_file(File.join(LocalPaths::Results, 'AllSNVs', 'YeastApobec', cancer_type.to_s, 'cancer.txt'))
    short_context_counts = short_context_counts_in_file(File.join(LocalPaths::Results, 'AllSNVs', 'YeastApobec', cancer_type.to_s, 'cancer.txt'))
    cancer_infos = ["#{cancer_type} (Yeast APOBEC sample)",
                    'N/A',
                    Configuration::YeastApobec::RandomShuffleFold[cancer_type],
                    'N/A',
                    Configuration::YeastApobec::FittingFoldShuffle[cancer_type]]
    matrix << [*cancer_infos,
              context_counts.each_value.inject(0, &:+),
              *possible_contexts.map{|context| context_counts[context] }]
    short_matrix << [*cancer_infos,
                    short_context_counts.each_value.inject(0, &:+),
                    *possible_short_contexts.map{|context| short_context_counts[context] }]
  }

  context_counts = context_counts_in_file(File.join(LocalPaths::Results, 'AllSNVs', 'NikZainal', 'cancer.txt'))
  short_context_counts = short_context_counts_in_file(File.join(LocalPaths::Results, 'AllSNVs', 'NikZainal', 'cancer.txt'))
  cancer_infos = ['Breast (NikZainal samples)',
                  Configuration::NikZainal::RandomGenomeFold,
                  Configuration::NikZainal::RandomShuffleFold,
                  Configuration::NikZainal::FittingFoldGenome,
                  Configuration::NikZainal::FittingFoldShuffle,
                ]
  matrix << [*cancer_infos,
              context_counts.each_value.inject(0, &:+),
              *possible_contexts.map{|context| context_counts[context] }]
  short_matrix << [*cancer_infos,
                  short_context_counts.each_value.inject(0, &:+),
                  *possible_short_contexts.map{|context| short_context_counts[context] }]

  File.open(File.join(LocalPaths::Results, 'motif_statistics/all_SNVs_sample_statistics.tsv'), 'w'){|fw|
    print_matrix(matrix.transpose, stream: fw)
  }
  File.open(File.join(LocalPaths::Results, 'motif_statistics/all_SNVs_sample_statistics_wo_mutation_direction.tsv'), 'w'){|fw|
    print_matrix(short_matrix.transpose, stream: fw)
  }

  #########
  rates_matrix = []
  rates_matrix << matrix.first
  matrix.drop(1).each{|row|
    total_count = row[cancer_infos.size]
    rates_matrix << [*row.first(cancer_infos.size), total_count, *row.drop(cancer_infos.size + 1).map{|count| (count.to_f / total_count).round(5) } ]
  }
  File.open(File.join(LocalPaths::Results, 'motif_statistics/all_SNVs_sample_statistics_rates.tsv'), 'w'){|fw|
    print_matrix(rates_matrix.transpose, stream: fw)
  }

  #########
  short_matrix_rates = []
  short_matrix_rates << short_matrix.first
  short_matrix.drop(1).each{|row|
    total_count = row[cancer_infos.size]
    short_matrix_rates << [*row.first(cancer_infos.size), total_count, *row.drop(cancer_infos.size + 1).map{|count| (count.to_f / total_count).round(5) } ]
  }
  File.open(File.join(LocalPaths::Results, 'motif_statistics/all_SNVs_sample_statistics_rates_wo_mutation_direction.tsv'), 'w'){|fw|
    print_matrix(short_matrix_rates.transpose, stream: fw)
  }
end

task :default => [:sample_statistics_all_SNVs]
