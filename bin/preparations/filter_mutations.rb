$:.unshift File.absolute_path('../../lib', __dir__)
require 'set'
require 'sequence_with_snp'
require 'optparse'
require 'mutation_context'
require 'perfectosape/results'
require 'breast_cancer_snv'
require 'import_information'

def snv_in_set(line, set)
  snv_normalized_name = line.split(' ', 2)[0].split('_')[0]
  set.include?(snv_normalized_name)
end

requested_mutation_region_types = nil
requested_mutation_contexts = nil
requested_cancer_samples = nil
invert_context_request = nil
show_possible_mutation_region_types_and_exit = false
show_possible_cancer_samples_and_exit = false
OptionParser.new do |opts|
  opts.banner = "Usage:\n" +
                "   #{opts.program_name} <SNV infos> <site fold-changes file> [options]\n" +
                "   <site fold-changes> | #{opts.program_name} <SNV infos> [options]"

  opts.separator 'Options:'
  opts.on('--mutation-types POSS_TYPES', "Specify possible mutation region types", "(e.g. `Intronic,Promoter`)") {|value|
    requested_mutation_region_types = RegionType.from_string(value)
  }

  opts.on('--contexts POSS_CONTEXTS',  "Specify possible mutation contexts (e.g. `TCN,NCG`)",
                                        "`N` represents any-nucleotide placeholder. `N`-s are independent",
                                        "Context is a 3nt string, with center representing a mutation" ) {|value|
    requested_mutation_contexts = value.upcase.split(',').map{|context| MutationContext.from_string(context, should_raise: false) }
  }

  opts.on('--cancer-samples POSSIBLE_SAMPLES', 'Specify cancer samples to retain.',
                                                'Format is comma separated. Be careful not to include spaces)',
                                                '(e.g. PD4194a,PD4198a') {|value|
    requested_cancer_samples = value.split(',').map(&:to_sym).to_set
  }

  opts.on('--show-possible-mutation-types', 'Show possible mutation types and exit') {
    show_possible_mutation_region_types_and_exit = true
  }

  opts.on('--show-possible-cancer-samples', 'Show possible cancer samples and exit') {
    show_possible_cancer_samples_and_exit = true
  }

  opts.on('--invert-context-request', 'Filter only mutations NOT matching specified context') {
    invert_context_request = true
  }
end.parse!(ARGV)

if requested_mutation_contexts
  requested_mutation_contexts = (requested_mutation_contexts + requested_mutation_contexts.map(&:revcomp)).uniq
else
  requested_mutation_contexts = [MutationContext.new('N','N','N', should_raise: false)]
end

raise 'Specify SNV infos file'  unless mutations_markup_filename = ARGV[0] # './source_data/SNV_infos.txt'

if show_possible_mutation_region_types_and_exit
  puts BreastCancerSNV.each_in_file(mutations_markup_filename).map(&:mutation_region_types).flat_map(&:features).inject(Set.new, &:merge).to_a
  exit
end

if show_possible_cancer_samples_and_exit
  puts BreastCancerSNV.each_in_file(mutations_markup_filename).flat_map(&:sample_id).inject(Set.new, &:<<).to_a
  exit
end

##########################

if invert_context_request
  context_filter = ->(snv){
    requested_mutation_contexts.none?{|requested_context| requested_context.match?(snv.context_before_snv_plus_strand) }
  }
else
  context_filter = ->(snv){
    requested_mutation_contexts.any?{|requested_context| requested_context.match?(snv.context_before_snv_plus_strand) }
  }
end

##########

if requested_cancer_samples
  sample_filter = ->(snv) {
    requested_cancer_samples.include?(snv.sample_id)
  }
else
  sample_filter = ->(snv){ true }
end

##########

if requested_mutation_region_types
  mutation_region_type_filter = ->(snv) {
    requested_mutation_region_types.features.intersect?(snv.mutation_region_types.features)
  }
else
  mutation_region_type_filter = ->(snv){ true }
end

##########

snvs_to_choose = BreastCancerSNV
  .each_in_file(mutations_markup_filename)
  .select(&context_filter)
  .select(&sample_filter)
  .select(&mutation_region_type_filter)
  .map(&:variant_id)
  .force
  .to_set

######################

if $stdin.tty?
  raise 'Specify site infos'  unless site_fold_changes_filename = ARGV[1] # './source_data/sites_cancer.txt'
  File.open(site_fold_changes_filename) do |f|
    f.each_line do |line|
      puts line  if snv_in_set(line, snvs_to_choose) # comment-line has normalized name starting with `#` so it's not in set
    end
  end
else
  $stdin.each_line do |line|
    puts line  if snv_in_set(line, snvs_to_choose) # comment-line has normalized name starting with `#` so it's not in set
  end
end
