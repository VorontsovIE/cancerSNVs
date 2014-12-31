$:.unshift File.absolute_path('../../lib', __dir__)
require 'support'
require 'mutation_features'

fasta_output_folder = 'results/motif_sites_fitted'

Dir.mkdir(fasta_output_folder)  unless Dir.exist?(fasta_output_folder)
Dir.glob(File.join(fasta_output_folder, '*')).select{|fn| File.file?(fn) }.each{|fn| File.delete(fn) }

motifs = Hash.new{|h,k| h[k] = [] }

mut_types = File.readlines('./source_data/SNV_infos.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type|
  intronic_mutation?(mut_type) || promoter_mutation?(mut_type)
}

MutatatedSiteInfo.each_site('./source_data/sites_cancer_any.txt').each do |info|
  next  unless info.pvalue_1 < 0.0005
  next  unless regulatory_mutation_names.include?(info.normalized_snp_name)

  motifs[info.motif_name] << info.seq_1
  if motifs[info.motif_name].size == 1000
    File.open(File.join(fasta_output_folder, "#{info.motif_name}.txt"),'a') do |fw|
      fw.puts motifs[info.motif_name]
    end
    motifs[info.motif_name] = []
  end
end

motifs.each do |motif_name, seqs|
  File.open(File.join(fasta_output_folder, "#{motif_name}.txt"),'a') do |fw|
    fw.puts seqs
  end
end
motifs
