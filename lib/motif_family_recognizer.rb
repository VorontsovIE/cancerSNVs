require_relative 'tf_ontology'
require_relative 'uniprot_reader'

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

  def subfamilies_by_multiple_uniprot_ids(uniprot_ids)
    uniprot_ids.flat_map{|uniprot_id|
      subfamilies_by_uniprot_id(uniprot_id)
    }.uniq
  end
end
