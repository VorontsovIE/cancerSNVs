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

def load_context_types(mutations_filename, mutations_markup_filename)
  snps_splitted = File.readlines(mutations_filename).map{|el| el.chomp.split("\t")}
  cpg_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| cpg_mutation?(sequence) }
  tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| tpc_mutation?(sequence) }
  not_cpg_tpc_names = mutation_names_by_mutation_context(snps_splitted){|mut_name, sequence| !tpc_mutation?(sequence) && !cpg_mutation?(sequence) }
  any_context_names = mutation_names_by_mutation_context(snps_splitted){ true }

  mut_types = File.readlines(mutations_markup_filename).drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
  intronic_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| intronic_mutation?(mut_type) }
  promoter_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type| promoter_mutation?(mut_type) }
  regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type|
    intronic_mutation?(mut_type) || promoter_mutation?(mut_type)
  }

  { tpcpg: tpc_names & cpg_names & regulatory_mutation_names,
    cpg: cpg_names & regulatory_mutation_names,
    tpc: tpc_names & regulatory_mutation_names,
    not_cpg_tpc: not_cpg_tpc_names & regulatory_mutation_names#,
    # any_context: any_context_names & regulatory_mutation_names
  }
end

def each_regulatory_site(mutation_filename, regulatory_mutation_names, context_types)
  each_mutated_site_info(mutation_filename) do |mutated_site_info|
    motif_name = mutated_site_info.motif_name
    pvalue_1 = mutated_site_info.pvalue_1
    snp_name = mutated_site_info.normalized_snp_name
    next  unless pvalue_1 <= 0.0005
    next  unless regulatory_mutation_names.include?(snp_name)

    contexts = context_type(snp_name, context_types)
    yield mutated_site_info, motif_name, contexts, pvalue_1
  end
end



# range = from...to
create_histogram = ->{  Histogram.new(1e-7, 0.0005, 1.0){|pvalue| - Math.log2(pvalue) }  }

context_aware = true

motif_names = File.readlines('source_data/motif_names.txt').map(&:strip)

context_types = load_context_types('source_data/SNPs.txt', 'source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt')

mutated_site_infos_cancer_filename = 'source_data/subsets/cancer_SNPs_regulatory_is_site.txt'
mutated_site_infos_shuffle_filename = 'source_data/subsets/shuffle_SNPs_regulatory_is_site.txt'
# mutated_site_infos_cancer_filename = 'source_data/cancer_SNPs.txt'
# mutated_site_infos_shuffle_filename = 'source_data/shuffle_SNPs.txt'

if context_aware
  fitters = ContextAwareMotifHistogramFitter.new(motif_names, context_types.keys, create_histogram.call)
else
  fitters = MotifHistogramFitter.new(motif_names, context_types.keys, create_histogram.call)
end


mut_types = File.readlines('source_data/SUBSTITUTIONS_13Apr2012_snz_promoter_markup2.txt').drop(1).map{|el| data = el.split("\t"); [data[0], data[17]] };
regulatory_mutation_names = mutation_names_by_mutation_type(mut_types){|mut_name, mut_type|
  intronic_mutation?(mut_type) || promoter_mutation?(mut_type)
}


each_regulatory_site(mutated_site_infos_cancer_filename, regulatory_mutation_names, context_types) do |mutated_site_info, motif_name, contexts, pvalue_1|
  fitters.add_element(motif_name, contexts, pvalue_1)
end

$stderr.puts "Loaded #{fitters.goal_total}"

num_iteration = 0
$stderr.puts "Start shuffle reading"
loop do
  num_iteration += 1

  each_regulatory_site(mutated_site_infos_shuffle_filename, regulatory_mutation_names, context_types) do |mutated_site_info, motif_name, contexts, pvalue_1|
    fitters.fit_element(motif_name, contexts, pvalue_1) { puts mutated_site_info.line }
  end

  raise StopIteration  if fitters.current_total >= fitters.goal_total
  fitters.print_discrepancies

  $stderr.puts "Loaded #{fitters.current_total}"
  raise StopIteration # Now we run the only iteration
end
