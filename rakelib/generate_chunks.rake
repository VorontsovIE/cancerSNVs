require 'shellwords'

def number_length(num)
  num.to_s.length
end

def common_length_number(num, max)
  '0' * (number_length(max) - number_length(num)) + num.to_s
end

def each_suffix_common_length(num)
  return enum_for(:each_suffix_common_length, num)  unless block_given?
  (1..num).each{|i|
    yield common_length_number(i, num)
  }
end

# Put links to neccessary resources into folder, generate script to run perfectosape on all files in dir
def prepare_core_folder(core_folder)
  mkdir_p core_folder
  ln_f LocalPaths::PerfectosAPE, File.join(core_folder, 'ape.jar')
  ln_f LocalPaths::MotifCollection, File.join(core_folder, 'motif_collection')
  ln_f LocalPaths::MotifThresholds,  File.join(core_folder, 'motif_thresholds')

  # create script from scratch
  File.open(File.join(core_folder, 'run_perfectosape.sh'), 'w') do |fw|
    fw.puts '#!/bin/bash'
  
    fw.puts 'cd "$(dirname "$0")"'
    Dir.glob(File.join(core_folder, '*.txt')).each do |filename|
      basename = File.basename(filename, '.txt')
      sites_filename = './sites_${basename}.txt'.shellescape
      fw.puts ['java', '${MEMORY_LIMIT}',
              '-cp ape.jar', 'ru.autosome.perfectosape.SNPScan',
              './motif_collection', filename.shellescape,
              '--fold-change-cutoff 1', '--precalc ./motif_thresholds', 
              '${EXPAND_FLANKS}', '--compact',
              "> #{sites_filename}"].join(" ")
    end
  end

  chmod 0755, File.join(core_folder, 'run_perfectosape.sh')
end

# Split SNV infos file into chunks and put them into folders related to cores
def split_file(filename, chunks_folder)
  mkdir_p chunks_folder
  basename = File.basename(filename, '.txt')
    sh 'split', "--number=l/#{Configuration::NumberOfCores}",
                '--numeric-suffixes=1',
                '--suffix-length=#{ number_length(Configuration::NumberOfCores) }',
                '--additional-suffix=.txt',
                filename,
                File.join(chunks_folder, "#{basename}_chunk_")
    
    each_suffix_common_length(Configuration::NumberOfCores) do |suffix|
      folder = File.join(chunks_folder, "core_#{suffix}")
      mkdir_p folder
      mv File.join(chunks_folder, "#{basename}_chunk_#{suffix}.txt"), 
         File.join(chunks_folder, "core_#{suffix}", "#{basename}.txt")
    end
end

# Generate scripts for concatenating results of parallel invocations
def create_concatenation_script(output_file, variants)
  mkdir_p File.dirname(output_file)
  File.open(output_file, 'w') do |fw|
    fw.puts '#!/bin/bash'
    fw.puts 'cd "$(dirname "$0")"'
    variants.each do |variant|
      fw.puts "cat /dev/null > #{filename}"
      # header
      sample_core_index = common_length_number(1, Configuration::NumberOfCores)
      fw.puts "grep -P ^# ./core_#{sample_core_index}/sites_#{variant}.txt  >>  ./sites_#{variant}.txt"

      # body
      each_suffix_common_length do |suffix|
        fw.puts "grep --invert-match -P ^# ./core_${suffix}/sites_#{variant}.txt  >>  ./sites_#{variant}.txt"
      end
      each_suffix_common_length do |suffix|
        fw.puts "rm ./core_#{suffix}/sites_#{variant}.txt"
      end
    end
  end
  chmod 0755, File.join(LocalPaths::Secondary::Chunks, 'concatenate_results.sh')
end

def create_script_for_multiple_invocation(output_file)
  mkdir_p File.dirname(output_file)
  # Prepare file to run all sequence chunks in parallel
  File.open(output_file, 'w') do |fw|
    fw.puts '#!/bin/bash'
    fw.puts "export EXPAND_FLANKS=\"--expand-region #{Configuration::ExpandFlanksLength}\""
    fw.puts "export MEMORY_LIMIT=\"#{Configuration::MemoryLimitOption}\""

    each_suffix_common_length(Configuration::NumberOfCores) do |suffix|
      fw.puts "./core_${suffix}/run_perfectosape.sh &"
    end
    fw.puts 'wait'
    fw.puts './concatenate_results.sh'
  end
  chmod 0755, File.join(LocalPaths::Secondary::Chunks, 'run_perfectosape_multithread.sh')
end

def prepare_chunks_for_sites(input_folder, output_folder)
  mkdir_p output_folder
  variants = Dir.glob(File.join(input_folder, '*.txt')).map{|fn| File.basename(fn, '.txt') }

  variants.each do |variant|
    split_file( File.join(input_folder, "#{variant}.txt"), File.join(output_folder) )
  end

  each_suffix_common_length(Configuration::NumberOfCores) do |suffix|
    prepare_core_folder File.join(output_folder, "core_#{suffix}")
  end

  create_script_for_multiple_invocation File.join(output_folder, 'run_perfectosape_multithread.sh')
  create_concatenation_script( File.join(output_folder, 'concatenate_results.sh'), variants )
end

namespace :preparations do
  desc 'Split sequences into equal chunks in order to run chunks in parallel'
  task generate_chunks: ['generate_chunks:Alexandrov', 'generate_chunks:NikZainal'] do
    File.open(File.join(LocalPaths::Secondary::Chunks, 'run_all.sh'), 'w') do |fw|
      fw.puts  File.join(LocalPaths::Secondary::Chunks, 'Alexandrov', 'run_all.sh')
      fw.puts  File.join(LocalPaths::Secondary::Chunks, 'NikZainal', 'run_perfectosape_multithread.sh')
    end
  end

  namespace 'generate_chunks' do
    task 'Alexandrov' do
      AlexandrovWholeGenomeCancers.each do |cancer_type|
        prepare_chunks_for_sites(
          File.join(LocalPaths::Secondary::SNVs, 'Alexandrov', cancer_type.to_s),
          File.join(LocalPaths::Secondary::Chunks, 'Alexandrov', cancer_type.to_s)
        )
      end
      File.open(File.join(LocalPaths::Secondary::Chunks, 'Alexandrov', 'run_all.sh'), 'w') do |fw|
        AlexandrovWholeGenomeCancers.each do |cancer_type|
          fw.puts  File.join(LocalPaths::Secondary::Chunks, 'Alexandrov', cancer_type.to_s, 'run_perfectosape_multithread.sh')
        end
      end
    end

    task 'NikZainal' do
      prepare_chunks_for_sites(
        File.join(LocalPaths::Secondary::NikZainalSNVs),
        File.join(LocalPaths::Secondary::Chunks, 'NikZainal')
      )
    end
  end
end
