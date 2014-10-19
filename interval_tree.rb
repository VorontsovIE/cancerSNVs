$:.unshift File.absolute_path('lib', __dir__)
$:.unshift '/home/ilya/iogen/projects/cage_analysis/lib/intervals/' # File.absolute_path('lib', __dir__)
require 'genome_region'
require 'mutation_features'
require 'support'

# class IntervalTree
#   attr_reader :chromosomes
#   def initialize(interval_object_pairs)
#     @chromosomes = interval_object_pairs
#       .group_by{|interval, object| "#{interval.chromosome},#{interval.strand}" }
#       .map{|chr, interval_object_pairs_subset|
#         # p interval_object_pairs_subset
#         intervals = interval_object_pairs_subset
#         .each_with_index
#         .sort_by{|(interval, object), ind|
#           interval.pos_end
#         }
#         [chr,intervals]
#       }.to_h#.tap{|x| p x}
#   end

#   def find(chr, strand, pos)
#     (interval, obj), ind = @chromosomes["#{chr},#{strand}"].bsearch{|(interval, obj), ind| interval.pos_end >= pos } # ?
#     @chromosomes["#{chr},#{strand}"].drop(ind).take_while{|(interval, obj), ind| interval.pos_end}
#     interval.include_position?(pos) ? [interval, obj] : nil
#     # return []  unless ind
#     # # select tss-es in upstream and in downstream (because promoter region can be assymetric)
#     # tsses[ [ind - 1, 0].max, 2 ].map{|tss,ind| tss}
#   end
# end
#
# interval_tree = IntervalTree.new(interval_object_pairs.first(20))
# p interval_tree.find('CHR_HSCHR6_MHC_QBL_CTG1', :+, 29887830)


def context_type(name_snp, context_types)
  context_pair = context_types.detect{|context_type, context_type_nameset| context_type_nameset.include?(name_snp) }
  context_pair && context_pair.first
end

interval_object_pairs_by_chromosome_non_reduced = Hash.new { |hash, key| hash[key] = [] }
File.open('source_data/gene_tss.txt') do |f|
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


snps_splitted = File.readlines('source_data/SNPs.txt').map{|el| el.chomp.split("\t")}

cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }
########

mut_infos = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el|
  data = el.chomp.split("\t")
  [data[0], [data[2], data[3], data[17]]] # mut_name, chr, pos, mut_type
}.to_h
intronic_mutation_names = mutation_names_by_mutation_type(mut_infos){|mut_name, (chr, pos, mut_type)| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_infos){|mut_name, (chr, pos, mut_type)| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_infos){|mut_name, (chr, pos, mut_type)| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }
############

context_types = { cpg: cpg_names & regulatory_mutation_names,
                  tpc: tpc_names & regulatory_mutation_names,
                  not_cpg_tpc: not_cpg_tpc_names & regulatory_mutation_names,
                  any_context: any_context_names & regulatory_mutation_names }


$stderr.puts "substitutions loaded"

ensg_to_hgnc = File.readlines('source_data/mart_export_ensg_hgnc.txt').drop(1).map{|line|
  ensg, enst, hgnc, hgnc_id = line.chomp.split("\t")
  [ensg, hgnc]
}.to_h

$stderr.puts "ensg-hgnc conversion loaded"


disrupting_mutations = Hash.new{|hsh,k| hsh[k] = []}
mutated_sites = each_mutation_infos('source_data/cancer_SNPs.txt').select(&disrupted_and_in_set_checker(regulatory_mutation_names))
mutated_sites.each{|line, name_snp, motif_name, fold_change, pvalue_1, pvalue_2|
  disrupting_mutations[motif_name] << [name_snp.split("_")[0], fold_change, pvalue_1, pvalue_2]
}

# $stderr.puts disrupting_mutations['HIF1A_si'].inspect

$stderr.puts "cancer_SNPs loaded"

# motif_for_analysis = ['HIF1A_si', 'TFE3_f1', 'CEBPG_si', 'SP3_f1', 'MYC_f1', 'ENOA_si']

motif_for_analysis = File.readlines('source_data/motif_names.txt').map(&:strip)

puts "motif\tmut_name\tchr\tpos\tmut_type\tcontext\tensg\tensg_to_hgnc[ensg]\tfold_change\tpvalue_1\tpvalue_2"
motif_for_analysis.each do |motif|
  ensgs = []
  disrupting_mutations[motif].each do |mut_name, fold_change, pvalue_1, pvalue_2|
    chr, pos, mut_type = mut_infos[mut_name]
    context = context_type(mut_name, context_types)
    pos = pos.to_i
    interval_object_pairs_by_chromosome[chr.to_s].select do |interval, ensg|
      interval.chromosome.to_s == chr.to_s && interval.include_position?(pos)
    end.each do |interval, ensg|
      ensgs << ensg
      puts "#{motif}\t#{mut_name}\t#{chr}\t#{pos}\t#{mut_type}\t#{context}\t#{ensg}\t#{ensg_to_hgnc[ensg]}\t#{fold_change}\t#{pvalue_1}\t#{pvalue_2}"
    end
  end
  $stderr.puts motif
  gene_ids = ensgs.map{|ensg| ensg_to_hgnc[ensg] || ensg }#.uniq
  # puts  "#{motif}:\n#{gene_ids.join(',')}\n#{gene_ids.uniq.sort.join(',')}"
  # puts  "#{motif}:\n#{gene_ids.uniq.sort.join(',')}"
end
