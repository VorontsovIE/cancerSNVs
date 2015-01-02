$:.unshift File.absolute_path('../../../lib', __dir__)
$:.unshift '/home/ilya/iogen/projects/cage_analysis/lib/intervals/'
require 'genome_region'
require 'mutation_context'
require 'import_information'
require 'breast_cancer_snv'

snv_infos_filename = './source_data/SNV_infos.txt'
snv_sequences_filename = './results/intermediate/SNV_sequences.txt'
sites_filename = './source_data/sites_cancer.txt'
motif_names_filename = './source_data/motif_names.txt'
ensg_hgnc_conversion_filename = './source_data/mart_export_ensg_hgnc.txt'
gene_tss_filename = './source_data/gene_tss.txt'

interval_object_pairs_by_chromosome_non_reduced = Hash.new { |hash, key| hash[key] = [] }
File.open(gene_tss_filename) do |f|
  f.each_line.lazy.drop(1).each do |line|
    ensg, _enst, _transcript_start, _transcript_end, gene_start, gene_end, chr, strand, _tss = line.chomp.split("\t")
    gene_start = gene_start.to_i
    gene_end = gene_end.to_i
    if strand == '1'
      interval = GenomeRegion.new(chr, :+, gene_start - 2000, [gene_end, gene_start + 500].max)
    elsif strand == '-1'
      interval = GenomeRegion.new(chr, :-, [gene_start, gene_end - 500].min, gene_end + 2000)
    end
    interval_object_pairs_by_chromosome_non_reduced[chr.to_s] << [interval, ensg]  if interval
  end
end

interval_object_pairs_by_chromosome = {}
interval_object_pairs_by_chromosome_non_reduced.each do |chr, intervals|
  interval_object_pairs_by_chromosome[chr] = intervals.uniq
end

contexts = {
  tpc: MutationContext.new('T','C','N'),
  cpg: MutationContext.new('N','C','G')
}.map{|name, context|
  [name, [context, context.revcomp]]
}.to_h

snvs_by_context = contexts.map {|context_name, requested_mutation_contexts|
  snv_names = matching_context_mutation_names(snv_sequences_filename, requested_mutation_contexts)
  [context_name, snv_names]
}.to_h

########

snvs = BreastCancerSNV.each_substitution_in_file(snv_infos_filename).map{|snv| [snv.variant_id, snv] }.to_h
regulatory_mutation_names = BreastCancerSNV.each_substitution_in_file(snv_infos_filename).select(&:regulatory?).map(&:variant_id).to_set
############


$stderr.puts "substitutions loaded"

ensg_to_hgnc = File.readlines(ensg_hgnc_conversion_filename).drop(1).map{|line|
  ensg, enst, hgnc, hgnc_id = line.chomp.split("\t")
  [ensg, hgnc]
}.to_h

$stderr.puts "ensg-hgnc conversion loaded"


disrupting_mutations = MutatatedSiteInfo.each_site(sites_filename).select{|site|
  regulatory_mutation_names.include?(site.normalized_snp_name) &&
  site.site_before_substitution?(pvalue_cutoff: 0.0005) &&
  site.disrupted?(fold_change_cutoff: 5)
}.group_by(&:motif_name)

$stderr.puts "Sites that were around SNVs loaded"

motif_for_analysis = File.readlines(motif_names_filename).map(&:strip)

puts "motif\tmut_name\tchr\tpos\tmut_type\tcontext\tensg\tensg_to_hgnc[ensg]\tfold_change\tpvalue_1\tpvalue_2"
motif_for_analysis.each do |motif|
  ensgs = []
  disrupting_mutations[motif].each do |mutated_site_info|
    mut_name = mutated_site_info.normalized_snp_name
    snv = snvs[mut_name]
    chr = snv.chr
    pos = snv.position
    mut_type = snv.mutation_types_string
    context = snvs_by_context.select{|context_name, snv_names| snv_names.include?(mut_name) }.keys.join(',')
    pos = pos.to_i
    interval_object_pairs_by_chromosome[chr.to_s].select do |interval, ensg|
      interval.chromosome.to_s == chr.to_s && interval.include_position?(pos)
    end.each do |interval, ensg|
      ensgs << ensg
      puts [motif, mut_name, chr, pos, mut_type, context, ensg, ensg_to_hgnc[ensg],
            mutated_site_info.fold_change, mutated_site_info.pvalue_1, mutated_site_info.pvalue_2 ].join("\t")
    end
  end
  $stderr.puts motif
  gene_ids = ensgs.map{|ensg| ensg_to_hgnc[ensg] || ensg }#.uniq
  # puts  "#{motif}:\n#{gene_ids.join(',')}\n#{gene_ids.uniq.sort.join(',')}"
  # puts  "#{motif}:\n#{gene_ids.uniq.sort.join(',')}"
end
