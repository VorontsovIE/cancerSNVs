$:.unshift File.absolute_path('lib', __dir__)
require 'set'
require 'statistical_significance'
require 'mutation_features'
require 'support'
require 'histogram'

class Range
  def round_to_s(rate)
    from = self.begin.round(rate)
    to = self.end.round(rate)
    exclude_end? ? "#{from}...#{to}" : "#{from}..#{to}"
  end
end


def context_type(name_snp, context_types)
  context_types.select{|context_type, context_type_nameset| context_type_nameset.include?(name_snp) }
                .map{|context_type, context_type_nameset| context_type }
end

def print_histogram_discrepancies(histogram_1, histogram_2, msg = nil)
  unless histogram_1.elements_total == histogram_2.elements_total
    $stderr.puts(msg)  if msg
    $stderr.puts "#{histogram_1.elements_total}; #{histogram_2.elements_total}"
    bin_1_iterator = histogram_1.each_bin
    bin_2_iterator = histogram_2.each_bin
    bin_1_iterator.zip(bin_2_iterator).each do |(bin_range_1, bin_1_count), (bin_range_2, bin_2_count)|
      $stderr.puts [bin_range_1.round_to_s(3), bin_1_count, bin_2_count].join("\t")  if bin_1_count != bin_2_count
    end
  end
end

def print_discrepancies_different_contexts(motif_names, context_types, histograms_1, histograms_2)
  motif_names.each do |motif|
    context_types.each_key do |context_type|
      print_histogram_discrepancies(histograms_1[motif][context_type], histograms_2[motif][context_type], "#{motif},#{context_type}")
    end
  end
end

def add_pvalue_to_histogram(histogram, pvalue)
  val = invlog(pvalue)
  histogram.add_element(val)
  histogram.in_range?(val) ? 1 : 0
end

def fit_pvalue_to_histogram(histogram_target, histogram_goal, pvalue)
  val = invlog(pvalue)
  count_target = histogram_target.bin_count_for_value(val)
  count_goal = histogram_goal.bin_count_for_value(val)
  # if count_target && count_target < count_goal
  if count_target < count_goal
    result = add_pvalue_to_histogram(histogram_target, pvalue)
    yield
    result
  else
    0
  end
end

def invlog(pvalue)
  -Math.log2(pvalue)
end

# range = from...to
create_histogram = ->{ Histogram.new(invlog(0.0005), invlog(1e-7), 1.0) }

use_different_contexts = true

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

snps_splitted = File.readlines('source_data/SNPs.txt').map{|el| el.chomp.split("\t")}
cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }

mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) || promoter_mutation?(mut_type) }

# context_types = { cpg: cpg_names & regulatory_mutation_names,
#                   tpc: tpc_names & regulatory_mutation_names,
#                   not_cpg_tpc: not_cpg_tpc_names & regulatory_mutation_names#,
#                   # any_context: any_context_names & regulatory_mutation_names
#                 }


context_types = { tpcpg: tpc_names & cpg_names & regulatory_mutation_names,
                  cpg: cpg_names & regulatory_mutation_names,
                  tpc: tpc_names & regulatory_mutation_names,
                  not_cpg_tpc: not_cpg_tpc_names & regulatory_mutation_names#,
                  # any_context: any_context_names & regulatory_mutation_names
                }



mutated_site_infos_cancer_filename = 'source_data/subsets/cancer_SNPs_regulatory_is_site.txt'
mutated_site_infos_shuffle_filename = 'source_data/subsets/shuffle_SNPs_regulatory_is_site.txt'
# mutated_site_infos_cancer_filename = 'source_data/cancer_SNPs.txt'
# mutated_site_infos_shuffle_filename = 'source_data/shuffle_SNPs.txt'

if use_different_contexts
  cancer_histograms_in_context = {}
  shuffle_histograms_in_context = {}

  motif_names.each do |motif_name|
    cancer_histograms_in_context[motif_name] = {}
    shuffle_histograms_in_context[motif_name] = {}
    context_types.each_key do |context_type|
      cancer_histograms_in_context[motif_name][context_type] = create_histogram.call
      shuffle_histograms_in_context[motif_name][context_type] = create_histogram.call
    end
  end
else
  cancer_histograms = {}
  shuffle_histograms = {}

  motif_names.each do |motif_name|
    cancer_histograms[motif_name] = create_histogram.call
    shuffle_histograms[motif_name] = create_histogram.call
  end
end


total = 0
each_mutated_site_info(mutated_site_infos_cancer_filename) do |mutated_site_info|
  motif_name = mutated_site_info.motif_name
  pvalue_1 = mutated_site_info.pvalue_1
  snp_name = mutated_site_info.normalized_snp_name
  next  unless pvalue_1 <= 0.0005
  next  unless regulatory_mutation_names.include?(snp_name)

  if use_different_contexts
    contexts = context_type(snp_name, context_types)
    raise 'Unknown context'  if contexts.empty?
    contexts.each do |context|
      total += add_pvalue_to_histogram(cancer_histograms_in_context[motif_name][context], pvalue_1)
    end
  else
    total += add_pvalue_to_histogram(cancer_histograms[motif_name], pvalue_1)
  end
end

$stderr.puts "Loaded #{total}"


new_total = 0
num_iteration = 0
$stderr.puts "Start shuffle reading"
loop do
  num_iteration += 1
  each_mutated_site_info(mutated_site_infos_shuffle_filename) do |mutated_site_info|
    motif_name = mutated_site_info.motif_name
    pvalue_1 = mutated_site_info.pvalue_1
    snp_name = mutated_site_info.normalized_snp_name
    next  unless pvalue_1 <= 0.0005
    next  unless regulatory_mutation_names.include?(snp_name)

    if use_different_contexts
      contexts = context_type(snp_name, context_types)
      raise 'Unknown context'  if contexts.empty?
      contexts.each do |context|
        new_total += fit_pvalue_to_histogram( shuffle_histograms_in_context[motif_name][context],
                                              cancer_histograms_in_context[motif_name][context],
                                              pvalue_1 ) { puts mutated_site_info.line }
      end
    else
      new_total += fit_pvalue_to_histogram( shuffle_histograms[motif_name],
                                            cancer_histograms[motif_name],
                                            pvalue_1 ) { puts mutated_site_info.line }
    end

    raise StopIteration  if new_total >= total
  end



  if use_different_contexts
    motif_names.each do |motif|
      context_types.each_key do |context_type|
        print_histogram_discrepancies(shuffle_histograms_in_context[motif][context_type],
                                      cancer_histograms_in_context[motif][context_type],
                                      "\n#{motif},#{context_type}")
      end
    end
  else
    motif_names.each do |motif|
      print_histogram_discrepancies(shuffle_histograms[motif],
                                    cancer_histograms[motif],
                                    "\n#{motif}")
    end
  end

  $stderr.puts "End of file reached #{new_total}"
  raise StopIteration # Now we run the only iteration
end

$stderr.puts "Finished #{new_total}"
