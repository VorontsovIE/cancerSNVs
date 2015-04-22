require 'interval_notation'
$:.unshift File.absolute_path('../../../lib', __dir__)
require 'mutation_context'
require 'import_information'
require 'snv_info'
require 'perfectosape/results'

snv_infos_filename = './source_data/SNV_infos.txt'
snv_sequences_filename = './results/intermediate/SNV_sequences.txt'
sites_filename = './source_data/sites_cancer.txt'
motif_names_filename = './source_data/motif_names.txt'
ensg_hgnc_conversion_filename = './unnecessary/source_data_unnecessary/mart_export_ensg_hgnc.txt'
gene_tss_filename = './source_data/deprecated/gene_tss_hg19.txt'

interval_object_pairs_by_chromosome_non_reduced = Hash.new { |hash, key| hash[key] = [] }
File.open(gene_tss_filename) do |f|
  f.each_line.lazy.drop(1).each do |line|
    ensg, _enst, _transcript_start, _transcript_end, gene_start, gene_end, chromosome, strand, _tss = line.chomp.split("\t")
    gene_start = gene_start.to_i
    gene_end = gene_end.to_i
    if strand == '1'
      interval = IntervalNotation::Syntax::Long.closed_open(gene_start - 2000, [gene_end, gene_start + 500].max)
    elsif strand == '-1'
      interval = IntervalNotation::Syntax::Long.closed_open([gene_start, gene_end - 500].min, gene_end + 2000)
    end
    interval_object_pairs_by_chromosome_non_reduced[chromosome.to_s] << [interval, ensg]  if interval
  end
end

interval_object_pairs_by_chromosome = {}
interval_object_pairs_by_chromosome_non_reduced.each do |chromosome, intervals|
  interval_object_pairs_by_chromosome[chromosome] = intervals.uniq
end

contexts = {
  tpc: MutationContext.new('T','C','N', should_raise: false),
  cpg: MutationContext.new('N','C','G', should_raise: false)
}.map{|name, context|
  [name, [context, context.revcomp]]
}.to_h

snvs_by_context = contexts.map {|context_name, requested_mutation_contexts|
  snv_names = matching_context_mutation_names(snv_sequences_filename, requested_mutation_contexts)
  [context_name, snv_names]
}.to_h

########

snvs = SNVInfo.each_in_file(snv_infos_filename).map{|snv| [snv.variant_id, snv] }.to_h
regulatory_mutation_names = SNVInfo.each_in_file(snv_infos_filename).select(&:regulatory?).map(&:variant_id).to_set
############


$stderr.puts "substitutions loaded"

ensg_to_hgnc = File.readlines(ensg_hgnc_conversion_filename).drop(1).map{|line|
  ensg, enst, hgnc, hgnc_id = line.chomp.split("\t")
  [ensg, hgnc]
}.to_h

$stderr.puts "ensg-hgnc conversion loaded"


disrupting_mutations = PerfectosAPE::Result.each_in_file(sites_filename).select{|site|
  regulatory_mutation_names.include?(site.normalized_snv_name) &&
  site.site_before_substitution?(pvalue_cutoff: 0.0005) &&
  site.disrupted?(fold_change_cutoff: 5)
}.group_by(&:motif_name)

$stderr.puts "Sites that were around SNVs loaded"

motif_for_analysis = File.readlines(motif_names_filename).map(&:strip)

puts ['motif','mut_name','chr','pos','mut_type','context','ensg','ensg_to_hgnc[ensg]','fold_change','pvalue_1','pvalue_2'].join("\t")
motif_for_analysis.each do |motif|
  ensgs = []
  disrupting_mutations[motif].each do |mutated_site_info|
    mut_name = mutated_site_info.normalized_snv_name
    snv = snvs[mut_name]
    chromosome = snv.chromosome
    pos = snv.position
    mut_type = snv.mutation_region_types
    context = snvs_by_context.select{|context_name, snv_names| snv_names.include?(mut_name) }.keys.join(',')
    pos = pos.to_i
    interval_object_pairs_by_chromosome[chromosome.to_s].select do |interval, ensg|
      interval.include_position?(pos)
    end.each do |interval, ensg|
      ensgs << ensg
      puts [motif, mut_name, chromosome, pos, mut_type, context, ensg, ensg_to_hgnc[ensg],
            mutated_site_info.fold_change, mutated_site_info.pvalue_1, mutated_site_info.pvalue_2 ].join("\t")
    end
  end
  $stderr.puts motif
  gene_ids = ensgs.map{|ensg| ensg_to_hgnc[ensg] || ensg }#.uniq
  # puts  "#{motif}:\n#{gene_ids.join(',')}\n#{gene_ids.uniq.sort.join(',')}"
  # puts  "#{motif}:\n#{gene_ids.uniq.sort.join(',')}"
end
