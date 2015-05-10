# require_relative 'preparations/prepare_alexandrov_data'
# require_relative 'preparations/prepare_nik_zainal_data'

# namespace 'preparations' do
#   desc 'Convert SNVs to a common format. Convert format, filter regulatory SNVs, remove duplicates'
#   task 'extractSNVs' => ['extractSNVs:AlexandrovEtAl', 'extractSNVs:NikZainalEtAl']
# end

# mkdir_p LocalPaths::Secondary::Sequences
# mkdir_p LocalPaths::Secondary::Chunks

$:.unshift File.join(LocalPaths::Root, 'lib')
require 'snv_info'
require 'region_type'
require 'genome_markup'
require 'data_import/breast_cancer_snv'
require 'data_import/cancer_mutations_from_alexandrov_et_al'
require 'data_import/sample_info'
require 'data_import/cancer_mutations_loading'

def markup_and_filter_SNVInfos(snv_infos_stream, genome_markup)
  snv_infos_stream
    .map{|snv_info|
      snv_info.marked_up(genome_markup)
    }
    .select(&:regulatory?)
    .map(&:in_pyrimidine_context)
    .to_a
    .uniq{|snv_info| # remove duplicates (in each sample separately)
      [snv_info.sample_id, snv_info.snv_sequence] # Here all SNVs are already in pyrimidine context
    }
end

def markup_and_filter_SNVInfos_to_file(snv_infos_stream, genome_markup, output_file:)
  File.open(output_file, 'w') do |output_stream|
    output_stream.puts(SNVInfo::HEADER)
    markup_and_filter_SNVInfos(snv_infos_stream, GENOME_MARKUP).each{|snv_info|
      output_stream.puts(snv_info)
    }
  end
end


namespace 'preparations' do
  desc 'Split sequences into equal chunks in order to run chunks in parallel'
  task 'split_into_chunks' do
    sh({'NUMBER_OF_CORES' => '8'}, 'prepare_sequences_for_perfectosape_run.sh')
  end

  namespace 'extractSNVs' do

    desc 'Convert Nik-Zainal\'s mutations to a common format. Convert format, filter regulatory only SNVs, remove duplicates.'
    task :NikZainalEtAl => [LocalPaths::Secondary::NikZainalSNVs]
    file LocalPaths::Secondary::NikZainalSNVs => [LocalPaths::Secondary::NikZainalSNVsOriginal] do
      GENOME_MARKUP = GENOME_MARKUP_LOADER.load_markup
      snv_infos_stream = BreastCancerSNV \
        .each_in_file(LocalPaths::Secondary::NikZainalSNVsOriginal) \
        .map{|snv| snv.to_snv_info(GENOME_READER, flank_length: 50) }
      markup_and_filter_SNVInfos_to_file(snv_infos_stream, GENOME_MARKUP, output_file: LocalPaths::Secondary::NikZainalSNVs)
    end

    desc 'Convert Alexandrov\'s mutations to a common format. Convert format, filter regulatory only SNVs, remove duplicates.'
    task :Alexandrov do
      GENOME_MARKUP = GENOME_MARKUP_LOADER.load_markup
      mutations_by_cancer = ALEXANDROV_MUTATIONS_LOADER.load

      mutations_by_cancer.each do |cancer_type, mutations|
        folder = File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type)
        mkdir_p folder

        snv_infos_stream = mutations
          .select(&:snv?)
          .map{|mutation|
            mutation.to_snv_info(GENOME_READER,
              cancer_type: cancer_type,
              variant_id: "#{mutation.sample_id};#{mutation.chromosome}:#{mutation.position_start}/#{mutation.after_substitution}",
              flank_length: 50)
          }
        markup_and_filter_SNVInfos_to_file(snv_infos_stream, GENOME_MARKUP, output_file: File.join(folder, "#{cancer_type}.txt"))
      end

    end
  end
end
