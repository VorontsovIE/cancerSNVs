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


interval_object_pairs_by_chromosome = Hash.new { |hash, key| hash[key] = [] }
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
    interval_object_pairs_by_chromosome[chr.to_s] << [interval, ensg]  if interval
  end
end


snps_splitted = File.readlines('source_data/SNPs.txt').map{|el| el.chomp.split("\t")}

cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }
########
mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }
############




regulatory_mutations = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|line|
  data = line.chomp.split("\t")
  [data[0], [data[2], data[3], data[17]]]
}.select{|mut_name,(chr,pos,type)|
  type.split(',').map(&:strip).include?('Intronic') || type.split(',').map(&:strip).include?('Promoter')
}.to_h
# mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };

$stderr.puts "substitutions loaded"

ensg_to_hgnc = File.readlines('source_data/mart_export_ensg_hgnc.txt').drop(1).map{|line|
  ensg, enst, hgnc, hgnc_id = line.chomp.split("\t")
  [ensg, hgnc]
}.to_h

$stderr.puts "ensg-hgnc conversion loaded"


disrupted_mutations = Hash.new{|hsh,k| hsh[k] = []}
File.open('source_data/cancer_SNPs.txt') do |f|
  f.each_line.drop(1).each do |line|
    mut_name, motif =  line.chomp.split("\t").first(2)
    disrupted_mutations[motif] << mut_name  if regulatory_mutations.has_key?(mut_name.split("_")[0])
  end
end


$stderr.puts "cancer_SNPs loaded"

# motif_for_analysis = ['HIF1A_si', 'TFE3_f1', 'CEBPG_si', 'SP3_f1', 'MYC_f1', 'ENOA_si']

motif_for_analysis = File.readlines('source_data/motif_names.txt').map(&:strip)

motif_for_analysis.each do |motif|
  ensgs = []
  disrupted_mutations[motif].each do |mut_name|
    chr, pos, type = regulatory_mutations[mut_name]
    pos = pos.to_i
    interval_object_pairs_by_chromosome[chr.to_s].select do |interval, object|
      interval.chromosome.to_s == chr.to_s && interval.include_position?(pos)
    end.each do |interval, object|
      ensgs << object
    end
  end
  $stderr.puts motif
  gene_ids = ensgs.map{|ensg| ensg_to_hgnc[ensg] || ensg }#.uniq
  puts  "#{motif}:\n#{gene_ids.join(',')}\n#{gene_ids.uniq.sort.join(',')}"
end
