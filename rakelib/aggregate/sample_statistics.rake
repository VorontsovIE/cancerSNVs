require_relative '../../lib/calculate_contexts'

directory File.join(LocalPaths::Results, 'motif_statistics')
desc 'Collect sample statistics'
task :sample_statistics_regulatory_SNVs => File.join(LocalPaths::Results, 'motif_statistics') do
  header = ['Sample', 'Genomic control fold', 'Shuffle control fold', 'Genomic control fitting fold', 'Shuffle control fitting fold', 'Mutations total']
  matrix = [ [*header, *possible_contexts] ]
  short_matrix = [ [*header, *possible_short_contexts] ]
  Configuration.WholeGenomeCancers.each{|cancer_type|
    context_counts = context_counts_in_file(File.join(LocalPaths::Results, 'SNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt'))
    short_context_counts = short_context_counts_in_file(File.join(LocalPaths::Results, 'SNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt'))
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

  File.open(File.join(LocalPaths::Results, 'motif_statistics/regulatory_SNVs_sample_statistics.tsv'), 'w'){|fw|
    print_matrix(matrix.transpose, stream: fw)
  }
  File.open(File.join(LocalPaths::Results, 'motif_statistics/regulatory_SNVs_sample_statistics_wo_mutation_direction.tsv'), 'w'){|fw|
    print_matrix(short_matrix.transpose, stream: fw)
  }

  #########
  cancer_infos_size = 5
  rates_matrix = []
  rates_matrix << matrix.first
  matrix.drop(1).each{|row|
    total_count = row[cancer_infos_size]
    rates_matrix << [*row.first(cancer_infos_size), total_count, *row.drop(cancer_infos_size + 1).map{|count| (count.to_f / total_count).round(5) } ]
  }
  File.open(File.join(LocalPaths::Results, 'motif_statistics/regulatory_SNVs_sample_statistics_rates.tsv'), 'w'){|fw|
    print_matrix(rates_matrix.transpose, stream: fw)
  }

  #########
  short_matrix_rates = []
  short_matrix_rates << short_matrix.first
  short_matrix.drop(1).each{|row|
    total_count = row[cancer_infos_size]
    short_matrix_rates << [*row.first(cancer_infos_size), total_count, *row.drop(cancer_infos_size + 1).map{|count| (count.to_f / total_count).round(5) } ]
  }
  File.open(File.join(LocalPaths::Results, 'motif_statistics/regulatory_SNVs_sample_statistics_rates_wo_mutation_direction.tsv'), 'w'){|fw|
    print_matrix(short_matrix_rates.transpose, stream: fw)
  }
end

desc 'Calculate number of samples in each cancer'
task :calculate_number_of_samples  do
  sample_counts = Configuration.WholeGenomeCancers.map{|cancer_type|
    cancer_fn = File.join(LocalPaths::Results, 'SNVs', 'Alexandrov', cancer_type.to_s, 'cancer.txt')
    samples = SNVInfo.each_in_file(cancer_fn).map{|snv_info|
      snv_info.variant_id.split(';').first
    }.to_a.uniq
    [cancer_type, samples.size]
  }.to_h
  File.open(File.join(LocalPaths::Results, 'motif_statistics/numbers_of_samples.txt'), 'w') do |fw|
    sample_counts.each{|cancer_type, num_samples|
      fw.puts "#{cancer_type}\t#{num_samples}"
    }
    fw.puts
    total_samples = sample_counts.values.inject(&:+)
    fw.puts "Total\t#{total_samples}"
  end
end

task :default => [:sample_statistics_regulatory_SNVs]
