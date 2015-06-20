require 'set'
require_relative '../../lib/uniprot_reader'
require_relative '../../lib/tf_ontology'

def motif_uniprots(filename)
  File.readlines(filename).drop(1).map{|line|
    motif, gene, quality, weight, human_uniprot, mouse_uniprot, consensus = line.chomp.split("\t")
    [motif, human_uniprot]
  }.to_h
end

class MotifFamilyRecognizer
  def initialize(deepness, uniprots_filename, tf_classification_filename)
    @deepness = deepness
    @uniprots_filename = uniprots_filename
    @tf_classification_filename = tf_classification_filename
  end

  def uniprot_ids_by_ac(uniprot_ac)
    @uniprot_ids_by_ac ||= begin
      mapping = UniprotEntry.each_in_file(@uniprots_filename)
                            .group_by(&:uniprot_ac)
                            .map{|uniprot_ac, uniprots|
                              [uniprot_ac, uniprots.map(&:uniprot_id)]
                            }.to_h
      mapping.default_proc = ->(h,k){h[k] = [] }
      mapping
    end
    @uniprot_ids_by_ac[uniprot_ac]
  end

  def uniprot_ids_by_multiple_acs(uniprot_acs)
    uniprot_acs.flat_map{|uniprot_ac| uniprot_ids_by_ac(uniprot_ac) }.uniq
  end

  def tf_classification
    @tf_classification ||= TFClassification.from_file(@tf_classification_filename)
  end

  def subtree_groups
    @subtree_groups ||= tf_classification.tf_groups(@deepness)
  end

  # In most cases Uniprot refers the only leaf, but in some cases it refers several leafs in different subtrees.
  # So we return an array of subfamilies
  def subfamilies_by_uniprot_id(uniprot_id)
    @subtree_root_by_uniprot_id ||= begin
      result = Hash.new{|h,k| h[k] = [] }

      subtree_groups.each{|group_root, group_leafs|
        uniprot_acs = group_leafs.flat_map(&:uniprot_ACs)
        uniprot_ids_by_multiple_acs(uniprot_acs).each{|uniprot_id|
          result[uniprot_id] << group_root
        }
      }
      result
    end
    @subtree_root_by_uniprot_id[uniprot_id]
  end
end


def collect_different_sample_statistics_gluing_subfamilies(sample_files, motif_family_recognizer, uniprot_by_motif)
  motif_uniprot_ids_by_sample = sample_files.map{|header, filename|
    [header, File.readlines(filename).map(&:strip).map{|motif| uniprot_by_motif[motif] }.compact.to_set]
  }.to_h

  motif_subfamilies_by_sample = motif_uniprot_ids_by_sample.map{|header, uniprot_ids|
    motif_subfamilies = uniprot_ids.flat_map{|uniprot_id|
      motif_family_recognizer.subfamilies_by_uniprot_id(uniprot_id).map(&:to_s)
    }.to_set
    [header, motif_subfamilies]
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

directory 'results/motif_statistics/aggregated/'

desc 'Aggregate common motifs over samples'
task 'aggregate_common_motifs' => ['results/motif_statistics/aggregated/'] do
  uniprot_by_motif = motif_uniprots(LocalPaths::Secondary::GeneInfos)
  motif_family_recognizer = MotifFamilyRecognizer.new(3,'source_data/human_uniprot.txt',
                                                        'source_data/TFOntologies/TFClass_human.obo')
  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      files = sample_files('results/motif_statistics/common/', 'any', protected_or_subjected, characteristic)
      File.open("results/motif_statistics/aggregated/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv", 'w') {|fw|
        print_matrix(collect_different_sample_statistics(files), stream: fw)
      }
      File.open("results/motif_statistics/aggregated/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued.csv", 'w') {|fw|
        print_matrix(collect_different_sample_statistics_gluing_subfamilies(files, motif_family_recognizer, uniprot_by_motif), stream: fw)
      }
    end
  end
end

directory 'results/motif_statistics/aggregated_wo_fitting/'

desc 'Aggregate common motifs over samples (w/o fitting)'
task 'aggregate_common_motifs_wo_fitting' => ['results/motif_statistics/aggregated_wo_fitting/'] do
  uniprot_by_motif = motif_uniprots(LocalPaths::Secondary::GeneInfos)
  motif_family_recognizer = MotifFamilyRecognizer.new(3,'source_data/human_uniprot.txt',
                                                        'source_data/TFOntologies/TFClass_human.obo')
  [:protected, :subjected].each do |protected_or_subjected|
    ['disruption', 'emergence', 'substitution-in-core'].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      files = sample_files('results/motif_statistics/common_wo_fitting/', 'any', protected_or_subjected, characteristic)
      File.open("results/motif_statistics/aggregated_wo_fitting/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv", 'w') {|fw|
        print_matrix(collect_different_sample_statistics(files), stream: fw)
      }
      File.open("results/motif_statistics/aggregated_wo_fitting/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued.csv", 'w') {|fw|
        print_matrix(collect_different_sample_statistics_gluing_subfamilies(files, motif_family_recognizer, uniprot_by_motif), stream: fw)
      }
    end
  end
end
