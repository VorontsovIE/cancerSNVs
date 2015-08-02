require 'set'
require_relative '../../lib/motif_family_recognizer'

def load_motif_qualities(filename)
  File.readlines(filename).drop(1).map{|line|
    motif, gene, quality, weight, human_uniprot, mouse_uniprot, consensus = line.chomp.split("\t")
    [motif, quality]
  }.to_h
end

def collect_different_sample_statistics_gluing_subfamilies(sample_files, motif_family_recognizer)
  motif_subfamilies_by_sample = sample_files.map{|sample, filename|
    motifs = File.readlines(filename).map(&:strip)
    # don't use subfamilies_by_multiple_motifs because it removes duplicated inner nodes resulting from different motifs, which we want to count downstream
    families = motifs.flat_map{|motif| motif_family_recognizer.subfamilies_by_motif(motif) }
    [sample, families]
  }.to_h

  term_occurences_matrix(motif_subfamilies_by_sample)
end

def collect_different_sample_statistics(sample_files)
  motifs_by_sample = sample_files.map{|header, filename|
    [header, File.readlines(filename).map(&:strip)]
  }
  term_occurences_matrix(motifs_by_sample)
end

def print_matrix(matrix, stream:)
  matrix.each do |row|
    stream.puts row.join("\t")
  end
end

def term_occurences_matrix(terms_by_sample)
  term_union = terms_by_sample.map{|sample, terms_in_sample| terms_in_sample }.inject(Set.new, :|).sort_by(&:to_s)
  matrix = []
  matrix << ['Sample', *term_union]
  terms_by_sample.each{|sample, terms_in_sample|
    term_count = term_union.map{|term|
      terms_in_sample.count{|el| el == term }
    }
    matrix << [ sample, *term_count.map{|count| count.zero? ? nil : count } ]
  }
  matrix.transpose
end

def with_motif_info_rows(matrix, motif_family_recognizers, motif_qualities)
  matrix_header = matrix.first.drop(1) # without "Motif" column
  matrix_data_by_motif = matrix.drop(1).map{|row| [row.first, row.drop(1)] }.to_h
  motifs = matrix_data_by_motif.keys

  result = []
  result << ['Motif', 'Motif quality', 'Motif families (level 3)', 'Motif families (level 4)', *matrix_header]
  motifs.each{|motif|
    families_3 = motif_family_recognizers[3].subfamilies_by_motif(motif).map(&:to_s).join('; ')
    families_4 = motif_family_recognizers[4].subfamilies_by_motif(motif).map(&:to_s).join('; ')
    result << [motif, motif_qualities[motif], families_3, families_4, *matrix_data_by_motif[motif]]
  }
  result
end

def sample_files(folder_common_motifs, context, protected_or_subjected, characteristic, with_yeast: true)
  result = AlexandrovWholeGenomeCancers.map{|sample|
    [sample, File.join(folder_common_motifs, 'Alexandrov', sample.to_s,
                        context.to_s, protected_or_subjected.to_s, characteristic.to_s, 'compared_to_each.txt') ]
  }
  if with_yeast
    result += YeastApobecSamples.map{|sample|
      [sample, File.join(folder_common_motifs, 'YeastApobec', sample.to_s,
                          context.to_s, protected_or_subjected.to_s, characteristic.to_s, 'compared_to_each.txt') ]
    }
  end
  result.to_h
end

def fitted_non_fitted_occurence_state(in_fitted, in_nonfitted)
  if in_fitted
    in_nonfitted ? 0 : 1
  else
    in_nonfitted ? -1 : nil
  end
end


def fitted_non_fitted_occurence_matrix(motifs_fitted_by_sample, motifs_nonfitted_by_sample, motif_family_recognizers, motif_qualities)
  all_motifs = (motifs_fitted_by_sample.values + motifs_nonfitted_by_sample.values).inject(&:|).to_a.sort # values are sets
  qualities = all_motifs.map{|motif|  motif_qualities[motif]  }
  motifs_subfamilies_3 = all_motifs.map{|motif|
    motif_family_recognizers[3].subfamilies_by_motif(motif).map(&:to_s).join('; ')
  }
  motifs_subfamilies_4 = all_motifs.map{|motif|
    motif_family_recognizers[4].subfamilies_by_motif(motif).map(&:to_s).join('; ')
  }
  results = []
  results << ['Motif', *all_motifs]
  results << ['Motif quality', *qualities]
  results << ['Motif families (level 3)', *motifs_subfamilies_3]
  results << ['Motif families (level 4)', *motifs_subfamilies_4]
  Configuration.sample_with_context_paths(with_nik_zainal: false, with_yeast: false).each_key{|sample_name|
    motifs_fitted = motifs_fitted_by_sample[sample_name]
    motifs_nonfitted = motifs_nonfitted_by_sample[sample_name]
    motifs_occurences = all_motifs.map{|motif|
      fitted_non_fitted_occurence_state(motifs_fitted.include?(motif), motifs_nonfitted.include?(motif))
    }
    results << [sample_name, *motifs_occurences]
  }

  results.transpose
end

directory 'results/motif_statistics/aggregated/'

desc 'Aggregate common motifs over samples'
task :aggregate_common_motifs => ['results/motif_statistics/aggregated/'] do
  output_folder = 'results/motif_statistics/aggregated/'
  motif_qualities = load_motif_qualities(LocalPaths::Secondary::GeneInfos)

  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      files = sample_files('results/motif_statistics/common/', 'any', protected_or_subjected, characteristic, with_yeast: false)
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.tsv"), 'w') {|fw|
        matrix = collect_different_sample_statistics(files)
        matrix_augmented = with_motif_info_rows(matrix, MOTIF_FAMILY_RECOGNIZERS, motif_qualities)
        print_matrix(matrix_augmented, stream: fw)
      }
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued_level_3.tsv"), 'w') {|fw|
        matrix = collect_different_sample_statistics_gluing_subfamilies(files, MOTIF_FAMILY_RECOGNIZERS[3])
        print_matrix(matrix, stream: fw)
      }
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued_level_4.tsv"), 'w') {|fw|
        matrix = collect_different_sample_statistics_gluing_subfamilies(files, MOTIF_FAMILY_RECOGNIZERS[4])
        print_matrix(matrix, stream: fw)
      }
    end
  end
end

directory 'results/motif_statistics/aggregated_wo_fitting/'

desc 'Aggregate common motifs over samples (w/o fitting)'
task :aggregate_common_motifs_wo_fitting => ['results/motif_statistics/aggregated_wo_fitting/'] do
  output_folder = 'results/motif_statistics/aggregated_wo_fitting/'
  motif_qualities = load_motif_qualities(LocalPaths::Secondary::GeneInfos)

  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      files = sample_files('results/motif_statistics/common_wo_fitting/', 'any', protected_or_subjected, characteristic, with_yeast: false)
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.tsv"), 'w') {|fw|
        matrix = collect_different_sample_statistics(files)
        matrix_augmented = with_motif_info_rows(matrix, MOTIF_FAMILY_RECOGNIZERS, motif_qualities)
        print_matrix(matrix_augmented, stream: fw)
      }
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued_level_3.tsv"), 'w') {|fw|
        matrix = collect_different_sample_statistics_gluing_subfamilies(files, MOTIF_FAMILY_RECOGNIZERS[3])
        print_matrix(matrix, stream: fw)
      }
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued_level_4.tsv"), 'w') {|fw|
        matrix = collect_different_sample_statistics_gluing_subfamilies(files, MOTIF_FAMILY_RECOGNIZERS[4])
        print_matrix(matrix, stream: fw)
      }
    end
  end
end

directory 'results/motif_statistics/aggregated_comparison'
desc 'Compare motif sets for experiment with and without fitting'
task :compare_fitted_to_unfitted => 'results/motif_statistics/aggregated_comparison' do
  output_folder = 'results/motif_statistics/aggregated_comparison/'
  motif_qualities = load_motif_qualities(LocalPaths::Secondary::GeneInfos)

  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      filename_last_part = File.join(protected_or_subjected.to_s, characteristic.to_s, 'compared_to_each.txt')

      motifs_fitted_by_sample = Configuration.sample_with_context_paths(with_nik_zainal: false, with_yeast: false).map{|sample_name, sample_path|
        motifs = File.readlines(
          File.join('results/motif_statistics/common/', sample_path, filename_last_part)
        ).map(&:chomp).to_set
        [sample_name, motifs]
      }.to_h

      motifs_nonfitted_by_sample = Configuration.sample_with_context_paths(with_nik_zainal: false, with_yeast: false).map{|sample_name, sample_path|
        motifs = File.readlines(
          File.join('results/motif_statistics/common_wo_fitting/', sample_path, filename_last_part)
        ).map(&:chomp).to_set
        [sample_name, motifs]
      }.to_h

      occurence_matrix = fitted_non_fitted_occurence_matrix(motifs_fitted_by_sample, motifs_nonfitted_by_sample, MOTIF_FAMILY_RECOGNIZERS, motif_qualities)

      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.tsv"), 'w') do |fw|
        print_matrix(occurence_matrix, stream: fw)
      end
    end
  end
end

task :default => [:aggregate_common_motifs, :aggregate_common_motifs_wo_fitting, :compare_fitted_to_unfitted]
