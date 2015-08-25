UniprotInfo = Struct.new(:uniprot_ac, :uniprot_id,  :primary_gene_name, :synonym_gene_names, :protein_names) do
  def all_gene_names
    [primary_gene_name] + synonym_gene_names
  end

  def self.from_string(line)
    uniprot_ac, uniprot_id,  primary_gene_name, synonym_gene_names, protein_names = line.chomp.split("\t")
    self.new(uniprot_ac, uniprot_id,  primary_gene_name, synonym_gene_names.split, protein_names)
  end

  def self.each_in_file(filename, &block)
    File.readlines(filename).drop(1).map{|line| self.from_string(line) }.each(&block)
  end
end

def read_uniprot_ids_by_motif(uniprot_by_motif_fn)
  File.readlines(uniprot_by_motif_fn).drop(1).map{|line|
    motif, human_uniprots = line.split("\t").first(2)
    [motif, human_uniprots.split(',')]
  }.to_h
end

def read_uniprot_info_by_id(uniprot_infos_fn)
  UniprotInfo.each_in_file(uniprot_infos_fn).to_a \
    .group_by(&:uniprot_id) \
    .map{|uniprot_id, infos|
      [uniprot_id, *infos]
    }.to_h
end

def read_uniprot_infos_by_motif(uniprot_by_motif_fn, uniprot_infos_fn)
  uniprot_ids_by_motif = read_uniprot_ids_by_motif(uniprot_by_motif_fn)
  uniprot_infos_by_id = read_uniprot_info_by_id(uniprot_infos_fn)
  uniprot_ids_by_motif.map{|motif, uniprot_ids|
    uniprot_infos = uniprot_ids.map{|uniprot_id|  uniprot_infos_by_id[uniprot_id]  }
    [motif, uniprot_infos]
  }.to_h
end
