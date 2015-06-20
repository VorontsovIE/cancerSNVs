require 'set'
require_relative '../../lib/motif_family_recognizer'

def motif_uniprots(filename)
  File.readlines(filename).drop(1).map{|line|
    motif, gene, quality, weight, human_uniprot, mouse_uniprot, consensus = line.chomp.split("\t")
    [motif, human_uniprot.split(',')]
  }.to_h
end

def load_motif_qualities(filename)
  File.readlines(filename).drop(1).map{|line|
    motif, gene, quality, weight, human_uniprot, mouse_uniprot, consensus = line.chomp.split("\t")
    [motif, quality]
  }.to_h
end

def motif_uniprots_in_file(filename, uniprots_by_motif)
  File.readlines(filename).map(&:strip).flat_map{|motif|
    uniprots_by_motif[motif]
  }.compact.to_set
end

def collect_different_sample_statistics_gluing_subfamilies(sample_files, motif_family_recognizer, uniprots_by_motif)
  motif_uniprot_ids_by_sample = sample_files.map{|sample, filename|
    [sample, motif_uniprots_in_file(filename, uniprots_by_motif)]
  }.to_h

  motif_subfamilies_by_sample = motif_uniprot_ids_by_sample.map{|sample, uniprot_ids|
    motif_subfamilies = motif_family_recognizer.subfamilies_by_multiple_uniprot_ids(uniprot_ids).map(&:to_s)
    [sample, motif_subfamilies]
  }.to_h

  term_occurences_matrix(motif_subfamilies_by_sample)
end


def collect_different_sample_statistics(sample_files)
  motifs_by_sample = sample_files.map{|header, filename|
    [header, File.readlines(filename).map(&:strip).to_set]
  }
  term_occurences_matrix(motifs_by_sample)
end

def print_matrix(matrix, stream:)
  matrix.each do |row|
    stream.puts row.join("\t")
  end
end

def term_occurences_matrix(terms_by_sample)
  term_union = terms_by_sample.map{|sample, terms| terms }.inject(Set.new, :|).sort
  matrix = []
  matrix << [nil, *term_union]
  terms_by_sample.each{|sample, terms|
    term_presence = term_union.map{|term| terms.include?(term) }
    matrix << [ sample, *term_presence.map{|present| present ? '1' : ''} ]
  }
  matrix.transpose
end

def with_motif_info_rows(matrix, motif_family_recognizer, uniprots_by_motif, motif_qualities)
  matrix_transposed = matrix.transpose
  motifs = matrix_transposed.first.drop(1)
  qualities = motifs.map{|motif|  motif_qualities[motif]  }
  motifs_subfamilies = motifs.map{|motif|
    uniprot_ids = uniprots_by_motif[motif]
    motif_family_recognizer.subfamilies_by_multiple_uniprot_ids(uniprot_ids).map(&:to_s).join(',')
  }

  [['Motif', *motifs], ['Motif quality', *qualities], ['Motif families', *motifs_subfamilies], *matrix_transposed.drop(1)].transpose
end

def sample_files(folder_common_motifs, context, protected_or_subjected, characteristic)
  ( AlexandrovWholeGenomeCancers.map{|sample|
    [sample, File.join(folder_common_motifs, 'Alexandrov', sample.to_s,
                        context.to_s, protected_or_subjected.to_s, characteristic.to_s, 'compared_to_each.txt') ]
  } +
  YeastApobecSamples.map{|sample|
    [sample, File.join(folder_common_motifs, 'YeastApobec', sample.to_s,
                        context.to_s, protected_or_subjected.to_s, characteristic.to_s, 'compared_to_each.txt') ]
  } ).to_h
end

def fitted_non_fitted_occurence_state(in_fitted, in_nonfitted)
  if in_fitted
    in_nonfitted ? 0 : 1
  else
    in_nonfitted ? -1 : nil
  end
end


def fitted_non_fitted_occurence_matrix(motifs_fitted_by_sample, motifs_nonfitted_by_sample)
  all_motifs = (motifs_fitted_by_sample.values + motifs_nonfitted_by_sample.values).inject(&:|).to_a.sort # values are sets
  results = []
  results << ['Motif', *all_motifs]
  Configuration.sample_paths.each_key{|sample_name|
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
task 'aggregate_common_motifs' => ['results/motif_statistics/aggregated/'] do
  output_folder = 'results/motif_statistics/aggregated/'
  uniprots_by_motif = motif_uniprots(LocalPaths::Secondary::GeneInfos)
  motif_qualities = load_motif_qualities(LocalPaths::Secondary::GeneInfos)
  motif_family_recognizer = MotifFamilyRecognizer.new(3,'source_data/human_uniprot.txt',
                                                        'source_data/TFOntologies/TFClass_human.obo')
  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      files = sample_files('results/motif_statistics/common/', 'any', protected_or_subjected, characteristic)
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv"), 'w') {|fw|
        matrix = collect_different_sample_statistics(files)
        matrix_augmented = with_motif_info_rows(matrix, motif_family_recognizer, uniprots_by_motif, motif_qualities)
        print_matrix(matrix_augmented, stream: fw)
      }
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued.csv"), 'w') {|fw|
        matrix = collect_different_sample_statistics_gluing_subfamilies(files, motif_family_recognizer, uniprots_by_motif)
        print_matrix(matrix, stream: fw)
      }
    end
  end
end

directory 'results/motif_statistics/aggregated_wo_fitting/'

desc 'Aggregate common motifs over samples (w/o fitting)'
task 'aggregate_common_motifs_wo_fitting' => ['results/motif_statistics/aggregated_wo_fitting/'] do
  output_folder = 'results/motif_statistics/aggregated_wo_fitting/'
  uniprots_by_motif = motif_uniprots(LocalPaths::Secondary::GeneInfos)
  motif_qualities = load_motif_qualities(LocalPaths::Secondary::GeneInfos)
  motif_family_recognizer = MotifFamilyRecognizer.new(3,'source_data/human_uniprot.txt',
                                                        'source_data/TFOntologies/TFClass_human.obo')
  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      files = sample_files('results/motif_statistics/common_wo_fitting/', 'any', protected_or_subjected, characteristic)
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv"), 'w') {|fw|
        matrix = collect_different_sample_statistics(files)
        matrix_augmented = with_motif_info_rows(matrix, motif_family_recognizer, uniprots_by_motif, motif_qualities)
        print_matrix(matrix_augmented, stream: fw)
      }
      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued.csv"), 'w') {|fw|
        matrix = collect_different_sample_statistics_gluing_subfamilies(files, motif_family_recognizer, uniprots_by_motif)
        print_matrix(matrix, stream: fw)
      }
    end
  end
end

directory 'results/motif_statistics/aggregated_comparison'

desc 'Compare motif sets for experiment with and without fitting'
task :compare_fitted_to_unfitted => 'results/motif_statistics/aggregated_comparison' do
  output_folder = 'results/motif_statistics/aggregated_comparison/'

  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      filename_last_part = File.join(protected_or_subjected.to_s, characteristic.to_s, 'compared_to_each.txt')

      motifs_fitted_by_sample = Configuration.sample_paths.map{|sample_name, sample_path|
        motifs = File.readlines(
          File.join('results/motif_statistics/common/', sample_path, filename_last_part)
        ).map(&:chomp).to_set
        [sample_name, motifs]
      }.to_h

      motifs_nonfitted_by_sample = Configuration.sample_paths.map{|sample_name, sample_path|
        motifs = File.readlines(
          File.join('results/motif_statistics/common_wo_fitting/', sample_path, filename_last_part)
        ).map(&:chomp).to_set
        [sample_name, motifs]
      }.to_h

      occurence_matrix = fitted_non_fitted_occurence_matrix(motifs_fitted_by_sample, motifs_nonfitted_by_sample)

      File.open(File.join(output_folder, "#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv"), 'w') do |fw|
        print_matrix(occurence_matrix, stream: fw)
      end
    end
  end

end
