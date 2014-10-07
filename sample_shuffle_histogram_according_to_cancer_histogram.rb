$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'histogram'


def invlog(pvalue)
  -Math.log2(pvalue)
end

def context_type(name_snp, context_types)
  context_pair = context_types.detect{|context_type, context_type_nameset| context_type_nameset.include?(name_snp) }
  context_pair && context_pair.first
end

create_histogram = ->{ Histogram.new(invlog(0.0005), invlog(1e-15), 1.0) }

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

# snps_splitted = File.readlines('source_data/SNPs.txt').map{|el| el.chomp.split("\t")}
# cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
# tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
# not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
# any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }

mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }

# context_types = { #cpg: cpg_names & regulatory_mutation_names,
#                   #tpc: tpc_names & regulatory_mutation_names,
#                   #not_cpg_tpc: not_cpg_tpc_names & regulatory_mutation_names,
#                   any_context: any_context_names & regulatory_mutation_names }

# cancer_histograms_in_context = Hash.new{|hsh,k| hsh[k] = {} }
# shuffle_histograms_in_context = Hash.new{|hsh,k| hsh[k] = {} }


# motif_names.each do |motif_name|
#   context_types.each do |context_type, context_type_nameset|
#     cancer_histograms_in_context[motif_name][context_type] = create_histogram.call
#     shuffle_histograms_in_context[motif_name][context_type] = create_histogram.call
#   end
# end

cancer_histograms = {}
shuffle_histograms = {}

motif_names.each do |motif_name|
  cancer_histograms[motif_name] = create_histogram.call
  shuffle_histograms[motif_name] = create_histogram.call
end

total = 0
each_mutation_infos('source_data/cancer_SNPs.txt') do |name_snp, motif_name, fold_change, pvalue_1, pvalue_2|
  next  unless pvalue_1 < 0.0005
  next  unless regulatory_mutation_names.include?(name_snp.split("_")[0])
  # context = context_type(name_snp.split("_")[0], context_types)
  # next  unless context
  # cancer_histograms_in_context[motif_name][ context ].add_element( invlog(pvalue_1) )
  cancer_histograms[motif_name].add_element( invlog(pvalue_1) )
  total += 1
end

$stderr.puts "Loaded #{total}"


new_total = 0
num_iteration = 0

$stderr.puts "Start shuffle reading"
loop do
  num_iteration += 1
  #!!!!!!!!!! HALF
  each_mutation_infos('source_data/half_shuffle_SNPs.txt') do |name_snp, motif_name, fold_change, pvalue_1, pvalue_2|
    next  unless pvalue_1 < 0.0005
    next  unless regulatory_mutation_names.include?(name_snp.split("_")[0])
    # context = context_type(name_snp.split("_")[0], context_types)
    # next  unless context
          
    # if shuffle_histograms_in_context[motif_name][context].elements_total < cancer_histograms_in_context[motif_name][context].elements_total
    #   val = invlog(pvalue_1)
    #   if shuffle_histograms_in_context[motif_name][context].bin_for_value(val) < cancer_histograms_in_context[motif_name][context].bin_for_value(val)
    #     shuffle_histograms_in_context[motif_name][context].add_element(val)
    #     new_total += 1
    #     # puts line
    #     # !!! add_to_output !!!
    #   end
    # end
    if shuffle_histograms[motif_name].elements_total < cancer_histograms[motif_name].elements_total
      val = invlog(pvalue_1)
      if shuffle_histograms[motif_name].bin_for_value(val) < cancer_histograms[motif_name].bin_for_value(val)
        shuffle_histograms[motif_name].add_element(val)
        new_total += 1
        # puts line
        # !!! add_to_output !!!
      end
    end

    raise StopIteration  if new_total >= total
  
  end
  $stderr.puts "End of file reached #{new_total}"
end

$stderr.puts "Finished #{new_total}"
