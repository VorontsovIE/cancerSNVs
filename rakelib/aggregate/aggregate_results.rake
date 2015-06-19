require 'set'
require_relative '../../lib/uniprot_reader'
require_relative '../../lib/tf_ontology'

def subfamilies_by_uniprot_ids(deepness, uniprots_filename, tf_classification_filename)
  result = Hash.new{|h,k| h[k] = [] }

  uniprots_by_AC = UniprotEntry.each_in_file(uniprots_filename).group_by(&:uniprot_ac)
  uniprots_by_AC.default_proc = ->(h,k){h[k] = [] }

  tf_classification = TFClassification.from_file(tf_classification_filename)
  subtree_groups = tf_classification.tf_groups(deepness)
  subtree_groups.each{|group_root, group_leafs|
    uniprot_acs = group_leafs.flat_map(&:uniprot_ACs)
    uniprot_ids = uniprot_acs.flat_map{|uniprot_ac| uniprots_by_AC[uniprot_ac] }.flat_map(&:uniprot_id).uniq
    uniprot_ids.each{|uniprot_id|
      result[uniprot_id] << group_root
    }
  }
  result
end

def motif_uniprots(filename)
  File.readlines(filename).drop(1).map{|line|
    motif, gene, quality, weight, human_uniprot, mouse_uniprot, consensus = line.chomp.split("\t")
    [motif, human_uniprot]
  }.to_h
end



def collect_different_sample_statistics_gluing_subfamilies(sample_files, stream:)
  uniprot_id_to_subtree_root = subfamilies_by_uniprot_ids(3, 'source_data/human_uniprot.txt',
                                                          'source_data/TFOntologies/TFClass_human.obo')
  uniprot_by_motif = motif_uniprots(LocalPaths::Secondary::GeneInfos)

  motif_uniprot_ids_by_sample = sample_files.map{|header, filename|
    [header, File.readlines(filename).map(&:strip).map{|motif| uniprot_by_motif[motif] }.to_set]
  }.to_h

  motif_uniprot_ids_by_sample.map{|header, uniprot_ids|
    motif_subfamilies = uniprot_ids.flat_map{|uniprot_id|
      uniprot_id_to_subtree_root[uniprot_id].map(&:to_s)
    }.to_set
    [header, motif_subfamilies]
  }.to_h

  print_term_occurences(motif_uniprot_ids_by_sample, stream: stream)
end


def collect_different_sample_statistics(sample_files, stream:)
  motifs_by_sample = sample_files.map{|header, filename|
    [header, File.readlines(filename).map(&:strip).to_set]
  }
  print_term_occurences(motifs_by_sample, stream: stream)
end

def print_term_occurences(terms_by_sample, stream:)
  term_union = terms_by_sample.map{|sample, terms| terms }.inject(Set.new, :|).sort
  stream.puts ['Sample', *term_union].join("\t")
  terms_by_sample.each{|sample, terms|
    term_presence = term_union.map{|term| terms.include?(term) }
    stream.puts [ sample, term_presence.map{|present| present ? 'X' : ''} ].join("\t")
  }
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
  [:protected, :subjected].each do |protected_or_subjected|
    [:disruption, :emergence, :substitution_in_core].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      File.open("results/motif_statistics/aggregated/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv", 'w') {|fw|
        collect_different_sample_statistics(sample_files('results/motif_statistics/common/', 'any', protected_or_subjected, characteristic), stream: fw)
      }
      File.open("results/motif_statistics/aggregated/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context_glued.csv", 'w') {|fw|
        collect_different_sample_statistics_gluing_subfamilies(sample_files('results/motif_statistics/common/', 'any', protected_or_subjected, characteristic), stream: fw)
      }
    end
  end
end

directory 'results/motif_statistics/aggregated_wo_fitting/'

desc 'Aggregate common motifs over samples (w/o fitting)'
task 'aggregate_common_motifs_wo_fitting' => ['results/motif_statistics/aggregated_wo_fitting/'] do
  [:protected, :subjected].each do |protected_or_subjected|
    [:disruption, :emergence, :substitution_in_core].each do |characteristic|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      File.open("results/motif_statistics/aggregated_wo_fitting/#{protected_or_subjected}_#{prep}_#{characteristic}_in_any_context.csv", 'w') {|fw|
        collect_different_sample_statistics(sample_files('results/motif_statistics/common_wo_fitting/', 'any', protected_or_subjected, characteristic), stream: fw)
      }
    end
  end
end
