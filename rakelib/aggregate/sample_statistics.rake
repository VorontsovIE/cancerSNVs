require_relative '../../lib/calculate_contexts'

desc 'Collect sample statistics'
task :sample_statistics do
  matrix = []
  matrix << ['Sample', 'Genome fitting fold', 'Random fitting fold', 'Mutations total', *possible_contexts]
  short_matrix = []
  short_matrix << ['Sample', 'Genome fitting fold', 'Random fitting fold', 'Mutations total', *possible_short_contexts]
  Configuration.getAlexandrovWholeGenomeCancers.each{|cancer_type|
    context_counts = context_counts_in_file(File.join('results/SNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt'))
    short_context_counts = short_context_counts_in_file(File.join('results/SNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt'))
    matrix << ["#{cancer_type} (Alexandrov et al. sample)", 
                Configuration::Alexandrov::FittingFoldGenome[cancer_type],
                Configuration::Alexandrov::FittingFoldShuffle[cancer_type],
                context_counts.each_value.inject(0, &:+),
                *possible_contexts.map{|context| context_counts[context] }]
    short_matrix << ["#{cancer_type} (Alexandrov et al. sample)",
                Configuration::Alexandrov::FittingFoldGenome[cancer_type],
                Configuration::Alexandrov::FittingFoldShuffle[cancer_type],
                short_context_counts.each_value.inject(0, &:+),
                *possible_short_contexts.map{|context| short_context_counts[context] }]
  }

  Configuration.getYeastApobecSamples.each{|cancer_type|
    context_counts = context_counts_in_file(File.join('results/SNVs', 'YeastApobec', cancer_type.to_s, 'cancer.txt'))
    short_context_counts = short_context_counts_in_file(File.join('results/SNVs', 'YeastApobec', cancer_type.to_s, 'cancer.txt'))
    matrix << ["#{cancer_type} (Yeast APOBEC sample)",
                'N/A',
                Configuration::YeastApobec::FittingFoldShuffle[cancer_type],
                context_counts.each_value.inject(0, &:+),
                *possible_contexts.map{|context| context_counts[context] }]
    short_matrix << ["#{cancer_type} (Yeast APOBEC sample)",
                'N/A',
                Configuration::YeastApobec::FittingFoldShuffle[cancer_type],
                short_context_counts.each_value.inject(0, &:+),
                *possible_short_contexts.map{|context| short_context_counts[context] }]
  }

  context_counts = context_counts_in_file(File.join('results/SNVs', 'NikZainal', 'cancer.txt'))
  short_context_counts = short_context_counts_in_file(File.join('results/SNVs', 'NikZainal', 'cancer.txt'))
  matrix << ['Breast (NikZainal samples)',
              Configuration::NikZainal::FittingFoldGenome,
              Configuration::NikZainal::FittingFoldShuffle,
              context_counts.each_value.inject(0, &:+),
              *possible_contexts.map{|context| context_counts[context] }]
  short_matrix << ['Breast (NikZainal samples)',
              Configuration::NikZainal::FittingFoldGenome,
              Configuration::NikZainal::FittingFoldShuffle,
              short_context_counts.each_value.inject(0, &:+),
              *possible_short_contexts.map{|context| short_context_counts[context] }]

  File.open('results/motif_statistics/sample_statistics.tsv', 'w'){|fw|
    print_matrix(matrix.transpose, stream: fw)
  }
  File.open('results/motif_statistics/sample_statistics_wo_mutation_direction.tsv', 'w'){|fw|
    print_matrix(short_matrix.transpose, stream: fw)
  }

  #########
  rates_matrix = []
  rates_matrix << matrix.first
  matrix.drop(1).each{|row|
    total_count = row[3]
    rates_matrix << [*row.first(3), total_count, *row.drop(4).map{|count| (count.to_f / total_count).round(5) } ]
  }
  File.open('results/motif_statistics/sample_statistics_rates.tsv', 'w'){|fw|
    print_matrix(rates_matrix.transpose, stream: fw)
  }

  #########
  short_matrix_rates = []
  short_matrix_rates << short_matrix.first
  short_matrix.drop(1).each{|row|
    total_count = row[3]
    short_matrix_rates << [*row.first(3), total_count, *row.drop(4).map{|count| (count.to_f / total_count).round(5) } ]
  }
  File.open('results/motif_statistics/sample_statistics_rates_wo_mutation_direction.tsv', 'w'){|fw|
    print_matrix(short_matrix_rates.transpose, stream: fw)
  }
end
