$:.unshift File.join(LocalPaths::Root, 'lib')
require 'snv_info'
require 'region_type'
require 'genome_markup'
require 'data_import/breast_cancer_snv'
require 'data_import/cancer_mutations_from_alexandrov_et_al'
require 'data_import/sample_info'
require 'data_import/cancer_mutations_loading'

def markup_dont_filter_SNVInfos(snv_infos_stream, genome_markup)
  snv_infos_stream
    .map{|snv_info|
      snv_info.marked_up(genome_markup)
    }
    .reject{|snv_info| snv_info.snv_sequence.to_s.match(/N/) }
    .to_a
    .uniq{|snv_info| # remove duplicates (in each sample separately)
      [snv_info.sample_id, snv_info.in_pyrimidine_context.snv_sequence]
    }
end

def markup_dont_filter_SNVInfos_to_file(snv_infos_stream, genome_markup, output_file:)
  mkdir_p File.dirname(output_file)
  File.open(output_file, 'w') do |output_stream|
    output_stream.puts(SNVInfo::HEADER)
    markup_dont_filter_SNVInfos(snv_infos_stream, genome_markup).each{|snv_info|
      output_stream.puts(snv_info)
    }
  end
end

# Nik-Zainal
folder = File.join(LocalPaths::Results, 'AllSNVs/NikZainal/')
directory folder
file File.join(LocalPaths::Results, 'AllSNVs/NikZainal/cancer.txt') => [LocalPaths::Secondary::NikZainalSNVsOriginal, folder] do
  snv_infos_stream = BreastCancerSNV \
    .each_in_file(LocalPaths::Secondary::NikZainalSNVsOriginal) \
    .map{|snv| snv.to_snv_info(GENOME_READER, flank_length: 50) }
  genome_markup = GENOME_MARKUP_LOADER.load_markup(dhs_accessible_filename: Configuration::DHS_BED_FILES[:Breast])
  markup_dont_filter_SNVInfos_to_file(snv_infos_stream, genome_markup, output_file: File.join(LocalPaths::Results, 'AllSNVs/NikZainal/cancer.txt'))
end

# Alexandrov
AlexandrovWholeGenomeCancers.each do |cancer_type|
  output_folder = File.join(LocalPaths::Results, 'AllSNVs', 'Alexandrov', cancer_type.to_s)
  cancer_filename = File.join(output_folder, 'cancer.txt')
  mutations_filename = File.join(LocalPaths::Secondary::Alexandrov::Mutations,
                                cancer_type.to_s,
                                "#{cancer_type}_clean_somatic_mutations_for_signature_analysis.txt")

  directory output_folder
  file  cancer_filename => [mutations_filename, LocalPaths::Secondary::Alexandrov::SamplesSummary, output_folder] do
    sample_infos = load_sample_infos(LocalPaths::Secondary::Alexandrov::SamplesSummary)
    snv_infos_stream = mutations_filered_by_sample(mutations_filename, sample_infos[cancer_type])
      .select(&:snv?)
      .map{|mutation|
        mutation.to_snv_info(GENOME_READER,
          cancer_type: cancer_type,
          variant_id: "#{mutation.sample_id};#{mutation.chromosome}:#{mutation.position_start}",
          flank_length: 50)
      }
    genome_markup = GENOME_MARKUP_LOADER.load_markup(dhs_accessible_filename: Configuration::DHS_BED_FILES[cancer_type])
    markup_dont_filter_SNVInfos_to_file(snv_infos_stream, genome_markup, output_file: cancer_filename)
  end
  task 'preparations:extractAllSNVs:Alexandrov' => cancer_filename
end

YeastApobecSamples.each do |sample|
  source_file = File.join('source_data/YeastApobec', "#{sample}.mfa")
  output_folder = File.join(LocalPaths::Results, 'AllSNVs/YeastApobec', sample.to_s)
  resulting_file = File.join(output_folder, 'cancer.txt')
  directory output_folder
  file resulting_file => [output_folder, source_file] do
    rm resulting_file  if File.exist?(resulting_file)
    cp source_file, resulting_file
  end
  task 'preparations:extractAllSNVs:YeastApobec' => resulting_file
end

namespace 'preparations' do
  desc 'Convert SNVs to proper place; DOESN\'T filter regulatory mutations as extractSNVs task do. But it also does remove duplicates. Necessary task to calculate number of non-regulatory mutations'
  task extractAllSNVs: ['extractAllSNVs:NikZainalEtAl', 'extractAllSNVs:Alexandrov', 'extractAllSNVs:YeastApobec']
  namespace 'extractAllSNVs' do
    desc 'Convert Nik-Zainal\'s mutations to a common format. Convert format, DOESN\'T filter regulatory only SNVs, remove duplicates.'
    task :NikZainalEtAl => [File.join(LocalPaths::Results, 'AllSNVs/NikZainal/cancer.txt')]

    desc 'Convert Alexandrov\'s mutations to a common format. Convert format, DOESN\'T filter regulatory only SNVs, remove duplicates.'
    task :Alexandrov

    desc 'Copy APOBEC mutations in yeast in proper format to necessary place.'
    task :YeastApobec
  end
end

def load_sample_infos(filename)
  $whole_genome_samples_by_cancer ||= whole_genome_samples_by_cancer(filename)
end
