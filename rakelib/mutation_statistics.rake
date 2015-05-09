# Statistics over cancer samples from Alexandrov et al.

def counts_by_property(mutations, property)
  mutations
    .group_by(&property)
    .map{|mutation_type, mutations_subset|
      [mutation_type, mutations_subset.size]
    }
    .sort_by{|mutation_type, size|
      size
    }.reverse.to_h
end

namespace :statistics do
  namespace :Alexandrov do
    desc 'Statistics of mutation types(SNV,DNV,TNV,indel) in different region types of different cancer types'
    task mutation_types: [LocalPaths::Results::Alexandrov::MutationTypes]
    file LocalPaths::Results::Alexandrov::MutationTypes do
      mutations_by_cancer_type = ALEXANDROV_MUTATIONS_LOADER.load
      File.open(LocalPaths::Results::Alexandrov::MutationTypesStatistics, 'w') do |fw|
        mutations_by_cancer_type.each do |cancer_type, mutations|
          fw.puts "> #{cancer_type}"

          RegionType.each_possible do |look_for_region_type|
            mutations_of_specified_region_type = mutations.select{|mutation|
              GENOME_MARKUP.get_region_type(mutation.chromosome, mutation.interval) == look_for_region_type
            }
            mutation_types_counts = counts_by_property(mutations_of_specified_region_type, :mutation_type)
            fw.puts [look_for_region_type.description, mutation_types_counts].join("\t")
          end
        end
      end
    end

    desc 'Statistics of contexts of SNVs in different cancer types of Alexandrov data set'
    task context_types: [LocalPaths::Results::Alexandrov::ContextStatistics]
    file LocalPaths::Results::Alexandrov::ContextStatistics do
      mutations_by_cancer_type = ALEXANDROV_MUTATIONS_LOADER.load
      File.open('./results/alexandrov_somatic_mutations_contexts.txt', 'w') do |fw|
        fw.puts '(promoter || intronic) && !kataegis'
        mutations_by_cancer_type.each do |cancer_type, mutations|
          fw.puts "> #{cancer_type}"

          regulatory_mutations = mutations.select(&:snv?).select{|mutation|
            GENOME_MARKUP.regulatory?(mutation.chromosome, mutation.interval)
          }

          contexts = regulatory_mutations.map{|mutation|
            mutation.to_snv_info(GENOME_READER, cancer_type: cancer_type)
          }.map(&:in_pyrimidine_context).map(&:context_before)

          fw.puts counts_by_property(context, &:itself)
        end
      end
    end
  end
end
