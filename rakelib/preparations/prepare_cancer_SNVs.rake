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
    .to_a
    .uniq{|snv_info| # remove duplicates (in each sample separately)
      [snv_info.sample_id, snv_info.in_pyrimidine_context.snv_sequence]
    }
end

def markup_and_filter_SNVInfos_to_file(snv_infos_stream, genome_markup, output_file:)
  mkdir_p File.dirname(output_file)
  File.open(output_file, 'w') do |output_stream|
    output_stream.puts(SNVInfo::HEADER)
    markup_and_filter_SNVInfos(snv_infos_stream, genome_markup).each{|snv_info|
      output_stream.puts(snv_info)
    }
  end
end

# Nik-Zainal
folder = File.join(LocalPaths::Secondary::SNVs, 'NikZainal')
directory folder
file LocalPaths::Secondary::NikZainalSNVs => [LocalPaths::Secondary::NikZainalSNVsOriginal, folder] do
  snv_infos_stream = BreastCancerSNV \
    .each_in_file(LocalPaths::Secondary::NikZainalSNVsOriginal) \
    .map{|snv| snv.to_snv_info(GENOME_READER, flank_length: 50) }
  markup_and_filter_SNVInfos_to_file(snv_infos_stream, GENOME_MARKUP_LOADER.load_markup, output_file: LocalPaths::Secondary::NikZainalSNVs)
end

# Alexandrov
AlexandrovWholeGenomeCancers.each do |cancer_type|
  folder = File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s)
  cancer_filename = File.join(folder, 'cancer.txt')
  mutations_filename = File.join(LocalPaths::Secondary::Alexandrov::Mutations,
                                cancer_type.to_s,
                                "#{cancer_type}_clean_somatic_mutations_for_signature_analysis.txt")

  directory folder
  file  cancer_filename => [mutations_filename, :load_sample_infos, folder] do
    snv_infos_stream = mutations_filered_by_sample(mutations_filename, WHOLE_GENOME_SAMPLES_BY_CANCER[cancer_type])
      .select(&:snv?)
      .map{|mutation|
        mutation.to_snv_info(GENOME_READER,
          cancer_type: cancer_type,
          variant_id: "#{mutation.sample_id};#{mutation.chromosome}:#{mutation.position_start}",
          flank_length: 50)
      }
    markup_and_filter_SNVInfos_to_file(snv_infos_stream, GENOME_MARKUP_LOADER.load_markup, output_file: cancer_filename)
  end
  task 'preparations:extractSNVs:Alexandrov' => cancer_filename
end

namespace 'preparations' do
  task extractSNVs: ['extractSNVs:NikZainalEtAl', 'extractSNVs:Alexandrov']
  namespace 'extractSNVs' do
    desc 'Convert Nik-Zainal\'s mutations to a common format. Convert format, filter regulatory only SNVs, remove duplicates.'
    task :NikZainalEtAl => [LocalPaths::Secondary::NikZainalSNVs]

    desc 'Convert Alexandrov\'s mutations to a common format. Convert format, filter regulatory only SNVs, remove duplicates.'
    task :Alexandrov
  end
end

task :load_sample_infos => LocalPaths::Secondary::Alexandrov::SamplesSummary do
  WHOLE_GENOME_SAMPLES_BY_CANCER ||= whole_genome_samples_by_cancer(LocalPaths::Secondary::Alexandrov::SamplesSummary)
end
