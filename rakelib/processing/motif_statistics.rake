def make_statistics_comparison_task(cancer_slices_folder:, random_slices_folder:, fitting_log:, output_file:, task_name:)
  output_folder = File.dirname(output_file)

  directory output_folder
  task task_name => "#{task_name}:full"
  task "#{task_name}:full" => output_file
  
  task output_file => [cancer_slices_folder, random_slices_folder, fitting_log, output_folder] do
    ruby 'summary.rb',  cancer_slices_folder,
                        random_slices_folder,
                        '--fitting-log', fitting_log,
                        '--correction', Configuration::CorrectionMethod,
                        '--expand-control-set', Configuration::ExpandControlSetFold.to_s,
                        {out: output_file}, {}
  end
end

def make_filtered_statistics_task(motif_statistics_file:, output_folder:, task_name:)
  task task_name => "#{task_name}:filtered"
  ['subjected', 'protected'].each do |subjected_or_protected|
    ['disruption',  'emergence', 'substitution-in-core'].each do |characteristic|
      output_subfolder = File.join(output_folder, subjected_or_protected.to_s, characteristic.to_s)
      directory output_subfolder
      output_filename = File.join(output_subfolder, File.basename(motif_statistics_file))
      task "#{task_name}:filtered" => output_filename
      file output_filename => [output_subfolder] do
        ruby 'filter_summary.rb', motif_statistics_file,
                                '--motif-qualities', 'A,B,C,D',
                                '--significance', 0.05.to_s,
                                '--characteristic', characteristic,
                                "--#{subjected_or_protected}",
                                {out: output_filename}, {}
      end
    end
  end
end

def make_common_motifs_tasks(folder:, random_datasets:, output_folder:, output_file:, task_name:)
  ['subjected', 'protected'].each do |subjected_or_protected|
    ['disruption',  'emergence', 'substitution-in-core'].each do |characteristic|
      input_files = random_datasets.map{|dataset| File.join(folder, subjected_or_protected.to_s, characteristic.to_s, "#{dataset}.csv") }
      output_subfolder = File.join(output_folder, subjected_or_protected.to_s, characteristic.to_s)
      output_file_fullname = File.join(output_subfolder, output_file)
      
      directory output_subfolder
      task task_name => output_file_fullname
      file output_file_fullname => [*input_files, output_subfolder] do
        ruby 'common_motifs.rb', *input_files, {out: output_file_fullname}, {}
      end
    end
  end
end


def make_all_common_motifs_tasks(folder:, output_folder:, configuration_module:, task_name:)
  task  task_name => ["#{task_name}:genome", "#{task_name}:shuffle", "#{task_name}:all"]

  make_common_motifs_tasks(
    folder: folder,
    random_datasets: configuration_module.const_get(:RandomGenomeDatasets),
    output_folder: output_folder,
    output_file: 'compared_to_each_genome.txt',
    task_name: "#{task_name}:genome"
  )

  make_common_motifs_tasks(
    folder: folder,
    random_datasets: configuration_module.const_get(:RandomShuffleDatasets),
    output_folder: output_folder,
    output_file: 'compared_to_each_shuffle.txt',
    task_name: "#{task_name}:shuffle"
  )


  make_common_motifs_tasks(
    folder: folder,
    random_datasets: configuration_module.const_get(:RandomDatasets),
    output_folder: output_folder,
    output_file: 'compared_to_each.txt',
    task_name: "#{task_name}:all"
  )
end

desc 'Calculate motif statistics (don\'t use in multitask mode)'
task 'motif_statistics' => ['motif_statistics:Alexandrov', 'motif_statistics:NikZainal', 'motif_statistics:YeastApobec']

task 'motif_statistics:Alexandrov'
AlexandrovWholeGenomeCancers.each do |cancer_type|
  task "motif_statistics:Alexandrov" => "motif_statistics:Alexandrov:#{cancer_type}"
  Configuration::Alexandrov.contexts_by_cancer_type(cancer_type).each do |context|
    task "motif_statistics:Alexandrov:#{cancer_type}" => "motif_statistics:Alexandrov:#{cancer_type}:#{context}"

    Configuration::Alexandrov::RandomDatasets.each do |random_dataset|
      task "motif_statistics:Alexandrov:#{cancer_type}:#{context}" => "motif_statistics:Alexandrov:#{cancer_type}:#{context}:#{random_dataset}"

      make_statistics_comparison_task(
        cancer_slices_folder: File.join(LocalPaths::Secondary::Slices, 'Alexandrov', cancer_type.to_s, context.to_s, 'cancer'),
        random_slices_folder: File.join(LocalPaths::Secondary::Slices, 'Alexandrov', cancer_type.to_s, context.to_s, random_dataset),
        fitting_log: File.join(LocalPaths::Secondary::LogFolder, 'Alexandrov', cancer_type.to_s, context.to_s, "#{random_dataset}.log"),
        output_file: File.join(LocalPaths::Secondary::MotifStatistics, 'full', 'Alexandrov', cancer_type.to_s, context.to_s, "#{random_dataset}.csv"),
        task_name: "motif_statistics:Alexandrov:#{cancer_type}:#{context}:#{random_dataset}"
      )

      make_filtered_statistics_task(motif_statistics_file: File.join(LocalPaths::Secondary::MotifStatistics, 'full', 'Alexandrov', cancer_type.to_s, context.to_s, "#{random_dataset}.csv"),
                                    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'Alexandrov', cancer_type.to_s, context.to_s),
                                    task_name: "motif_statistics:Alexandrov:#{cancer_type}:#{context}:#{random_dataset}")
    end

    task "motif_statistics:Alexandrov:#{cancer_type}:#{context}" => "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs"
    make_all_common_motifs_tasks(
      folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'Alexandrov', cancer_type.to_s, context.to_s),
      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'Alexandrov', cancer_type.to_s, context.to_s),
      configuration_module: Configuration::Alexandrov,
      task_name: "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs"
    )
  end
end

task 'motif_statistics:NikZainal'
Configuration::NikZainalContexts.each do |context|
  task "motif_statistics:NikZainal" => "motif_statistics:NikZainal:#{context}"

  Configuration::NikZainal::RandomDatasets.each do |random_dataset|
    task "motif_statistics:NikZainal:#{context}" => "motif_statistics:NikZainal:#{context}:#{random_dataset}"

    make_statistics_comparison_task(
      cancer_slices_folder: File.join(LocalPaths::Secondary::Slices, 'NikZainal', context.to_s, 'cancer'),
      random_slices_folder: File.join(LocalPaths::Secondary::Slices, 'NikZainal', context.to_s, random_dataset),
      fitting_log: File.join(LocalPaths::Secondary::LogFolder, 'NikZainal', context.to_s, "#{random_dataset}.log"),
      output_file: File.join(LocalPaths::Secondary::MotifStatistics, 'full', 'NikZainal', context.to_s, "#{random_dataset}.csv"),
      task_name: "motif_statistics:NikZainal:#{context}:#{random_dataset}"
    )

    make_filtered_statistics_task(motif_statistics_file: File.join(LocalPaths::Secondary::MotifStatistics, 'full', 'NikZainal', context.to_s, "#{random_dataset}.csv"),
                                  output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'NikZainal', context.to_s),
                                  task_name: "motif_statistics:NikZainal:#{context}:#{random_dataset}")
  end

  task "motif_statistics:NikZainal:#{context}" => "motif_statistics:NikZainal:#{context}:common_motifs"
  make_all_common_motifs_tasks(
    folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'NikZainal', context.to_s),
    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'NikZainal', context.to_s),
    configuration_module: Configuration::NikZainal,
    task_name: "motif_statistics:NikZainal:#{context}:common_motifs"
  )
end

task 'motif_statistics:YeastApobec'
YeastApobecSamples.each do |sample|
  task "motif_statistics:YeastApobec" => "motif_statistics:YeastApobec:#{sample}"
  Configuration::YeastApobec.contexts_by_cancer_type(sample).each do |context|
    task "motif_statistics:YeastApobec:#{sample}" => "motif_statistics:YeastApobec:#{sample}:#{context}"

    Configuration::YeastApobec::RandomDatasets.each do |random_dataset|
      task "motif_statistics:YeastApobec:#{sample}:#{context}" => "motif_statistics:YeastApobec:#{sample}:#{context}:#{random_dataset}"

      make_statistics_comparison_task(
        cancer_slices_folder: File.join(LocalPaths::Secondary::Slices, 'YeastApobec', sample.to_s, context.to_s, 'cancer'),
        random_slices_folder: File.join(LocalPaths::Secondary::Slices, 'YeastApobec', sample.to_s, context.to_s, random_dataset),
        fitting_log: File.join(LocalPaths::Secondary::LogFolder, 'YeastApobec', sample.to_s, context.to_s, "#{random_dataset}.log"),
        output_file: File.join(LocalPaths::Secondary::MotifStatistics, 'full', 'YeastApobec', sample.to_s, context.to_s, "#{random_dataset}.csv"),
        task_name: "motif_statistics:YeastApobec:#{sample}:#{context}:#{random_dataset}"
      )

      make_filtered_statistics_task(motif_statistics_file: File.join(LocalPaths::Secondary::MotifStatistics, 'full', 'YeastApobec', sample.to_s, context.to_s, "#{random_dataset}.csv"),
                                    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'YeastApobec', sample.to_s, context.to_s),
                                    task_name: "motif_statistics:YeastApobec:#{sample}:#{context}:#{random_dataset}")
    end

    task "motif_statistics:YeastApobec:#{sample}:#{context}" => "motif_statistics:YeastApobec:#{sample}:#{context}:common_motifs"
    make_all_common_motifs_tasks(
      folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'YeastApobec', sample.to_s, context.to_s),
      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'YeastApobec', sample.to_s, context.to_s),
      configuration_module: Configuration::YeastApobec,
      task_name: "motif_statistics:YeastApobec:#{sample}:#{context}:common_motifs"
    )
  end
end
