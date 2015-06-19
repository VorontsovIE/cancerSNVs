def make_statistics_comparison_task(cancer_slices_folder:, random_slices_folder:, fitting_log:, output_file:, task_name:)
  output_folder = File.dirname(output_file)

  directory output_folder
  task task_name => "#{task_name}:full"
  task "#{task_name}:full" => output_file
  fitting_log_option = fitting_log ? ['--fitting-log', fitting_log] : []
  task output_file => [cancer_slices_folder, random_slices_folder, fitting_log, output_folder].compact do
    ruby 'summary.rb',  cancer_slices_folder,
                        random_slices_folder,
                        *fitting_log_option,
                        '--correction', Configuration::CorrectionMethod,
                        '--expand-control-set', Configuration::ExpandControlSetFold.to_s,
                        {out: output_file}, {}
  end
end

def make_filtered_statistics_task(motif_statistics_file:, output_folder:, task_name:)
  task task_name do
    ['subjected', 'protected'].each do |subjected_or_protected|
      ['disruption',  'emergence', 'substitution-in-core'].each do |characteristic|
        output_subfolder = File.join(output_folder, subjected_or_protected.to_s, characteristic.to_s)
        mkdir_p output_subfolder
        output_filename = File.join(output_subfolder, File.basename(motif_statistics_file))
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
  task task_name do
    ['subjected', 'protected'].each do |subjected_or_protected|
      ['disruption',  'emergence', 'substitution-in-core'].each do |characteristic|
        input_files = random_datasets.map{|dataset| File.join(folder, subjected_or_protected.to_s, characteristic.to_s, "#{dataset}.csv") }
        output_subfolder = File.join(output_folder, subjected_or_protected.to_s, characteristic.to_s)
        output_file_fullname = File.join(output_subfolder, output_file)

        mkdir_p output_subfolder
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

desc 'Calculate motif statistics; 1-st stage'
task 'motif_statistics' => ['motif_statistics:Alexandrov', 'motif_statistics:NikZainal', 'motif_statistics:YeastApobec']

desc 'Calculate motif statistics (w/o fitting); 1-st stage'
task 'motif_statistics_wo_fitting' => ['motif_statistics_wo_fitting:Alexandrov', 'motif_statistics_wo_fitting:NikZainal', 'motif_statistics_wo_fitting:YeastApobec']

desc 'Process motif statistics; 2-nd stage'
task 'filtered_motif_statistics' => ['filtered_motif_statistics:Alexandrov', 'filtered_motif_statistics:NikZainal', 'filtered_motif_statistics:YeastApobec']

desc 'Process motif statistics (w/o fitting); 2-nd stage'
task 'filtered_motif_statistics_wo_fitting' => ['filtered_motif_statistics_wo_fitting:Alexandrov', 'filtered_motif_statistics_wo_fitting:NikZainal', 'filtered_motif_statistics_wo_fitting:YeastApobec']

desc 'Collect common motifs from motif statistics; 3-rd stage'
task 'common_motif_statistics' => ['common_motif_statistics:Alexandrov', 'common_motif_statistics:NikZainal', 'common_motif_statistics:YeastApobec']

desc 'Collect common motifs from motif statistics (w/o fitting); 3-rd stage'
task 'common_motif_statistics_wo_fitting' => ['common_motif_statistics_wo_fitting:Alexandrov', 'common_motif_statistics_wo_fitting:NikZainal', 'common_motif_statistics_wo_fitting:YeastApobec']

def prefixed_motif_statistics_task(task, dependencies, perform_calculation: true, filtering: true, common_motifs: true)
  task_prefixes = []
  task_prefixes += ['motif_statistics', 'motif_statistics_wo_fitting']  if perform_calculation
  task_prefixes += ['filtered_motif_statistics', 'filtered_motif_statistics_wo_fitting']  if filtering
  task_prefixes += ['common_motif_statistics', 'common_motif_statistics_wo_fitting']  if common_motifs

  task_prefixes.each{|task_prefix|
    task "#{task_prefix}:#{task}" => Array(dependencies).map{|dependency| "#{task_prefix}:#{dependency}" }
  }
end

fitting_wo_fitting_settings = [
  ['motif_statistics', LocalPaths::Secondary::Slices, LocalPaths::Secondary::LogFolder, 'full', 'filtered', 'common'],
  ['motif_statistics_wo_fitting', 'results/motif_statistics/slices_wo_fitting', nil, 'full_wo_fitting', 'filtered_wo_fitting', 'common_wo_fitting']
]

fitting_wo_fitting_settings.each do |task_prefix, slices_folder, log_folder, full_folder, filtered_folder, common_folder|
  ####################################
  prefixed_motif_statistics_task 'Alexandrov', []
  AlexandrovWholeGenomeCancers.each do |cancer_type|
    prefixed_motif_statistics_task 'Alexandrov', ["Alexandrov:#{cancer_type}"]
    Configuration::Alexandrov.contexts_by_cancer_type(cancer_type).each do |context|
      prefixed_motif_statistics_task "Alexandrov:#{cancer_type}", ["Alexandrov:#{cancer_type}:#{context}"]

      Configuration::Alexandrov::RandomDatasets.each do |random_dataset|
        prefixed_motif_statistics_task "Alexandrov:#{cancer_type}:#{context}", "Alexandrov:#{cancer_type}:#{context}:#{random_dataset}", common_motifs: false
        make_statistics_comparison_task(
          cancer_slices_folder: File.join(slices_folder, 'Alexandrov', cancer_type.to_s, context.to_s, 'cancer'),
          random_slices_folder: File.join(slices_folder, 'Alexandrov', cancer_type.to_s, context.to_s, random_dataset),
          fitting_log: log_folder && File.join(log_folder, 'Alexandrov', cancer_type.to_s, context.to_s, "#{random_dataset}.log"),
          output_file: File.join(LocalPaths::Secondary::MotifStatistics, full_folder, 'Alexandrov', cancer_type.to_s, context.to_s, "#{random_dataset}.csv"),
          task_name: "#{task_prefix}:Alexandrov:#{cancer_type}:#{context}:#{random_dataset}"
        )

        make_filtered_statistics_task(motif_statistics_file: File.join(LocalPaths::Secondary::MotifStatistics, full_folder, 'Alexandrov', cancer_type.to_s, context.to_s, "#{random_dataset}.csv"),
                                      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, filtered_folder, 'Alexandrov', cancer_type.to_s, context.to_s),
                                      task_name: "filtered_#{task_prefix}:Alexandrov:#{cancer_type}:#{context}:#{random_dataset}")
      end

      make_all_common_motifs_tasks(
        folder: File.join(LocalPaths::Secondary::MotifStatistics, filtered_folder, 'Alexandrov', cancer_type.to_s, context.to_s),
        output_folder: File.join(LocalPaths::Secondary::MotifStatistics, common_folder, 'Alexandrov', cancer_type.to_s, context.to_s),
        configuration_module: Configuration::Alexandrov,
        task_name: "common_#{task_prefix}:Alexandrov:#{cancer_type}:#{context}"
      )
    end
  end

  ####################################

  prefixed_motif_statistics_task 'NikZainal', []
  Configuration::NikZainalContexts.each do |context|
    prefixed_motif_statistics_task "NikZainal", ["NikZainal:#{context}"]
    Configuration::NikZainal::RandomDatasets.each do |random_dataset|
      prefixed_motif_statistics_task "NikZainal:#{context}", ["NikZainal:#{context}:#{random_dataset}"], common_motifs: false

      make_statistics_comparison_task(
        cancer_slices_folder: File.join(slices_folder, 'NikZainal', context.to_s, 'cancer'),
        random_slices_folder: File.join(slices_folder, 'NikZainal', context.to_s, random_dataset),
        fitting_log: log_folder && File.join(log_folder, 'NikZainal', context.to_s, "#{random_dataset}.log"),
        output_file: File.join(LocalPaths::Secondary::MotifStatistics, full_folder, 'NikZainal', context.to_s, "#{random_dataset}.csv"),
        task_name: "#{task_prefix}:NikZainal:#{context}:#{random_dataset}"
      )

      make_filtered_statistics_task(motif_statistics_file: File.join(LocalPaths::Secondary::MotifStatistics, full_folder, 'NikZainal', context.to_s, "#{random_dataset}.csv"),
                                    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, filtered_folder, 'NikZainal', context.to_s),
                                    task_name: "filtered_#{task_prefix}:NikZainal:#{context}:#{random_dataset}")
    end

    make_all_common_motifs_tasks(
      folder: File.join(LocalPaths::Secondary::MotifStatistics, filtered_folder, 'NikZainal', context.to_s),
      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, common_folder, 'NikZainal', context.to_s),
      configuration_module: Configuration::NikZainal,
      task_name: "common_#{task_prefix}:NikZainal:#{context}"
    )
  end

  ####################################

  prefixed_motif_statistics_task 'YeastApobec', []
  YeastApobecSamples.each do |sample|
    prefixed_motif_statistics_task "YeastApobec", ["YeastApobec:#{sample}"]
    Configuration::YeastApobec.contexts_by_cancer_type(sample).each do |context|
      prefixed_motif_statistics_task "YeastApobec:#{sample}", ["YeastApobec:#{sample}:#{context}"]
      Configuration::YeastApobec::RandomDatasets.each do |random_dataset|
        prefixed_motif_statistics_task "YeastApobec:#{sample}:#{context}", ["YeastApobec:#{sample}:#{context}:#{random_dataset}"], common_motifs: false

        make_statistics_comparison_task(
          cancer_slices_folder: File.join(slices_folder, 'YeastApobec', sample.to_s, context.to_s, 'cancer'),
          random_slices_folder: File.join(slices_folder, 'YeastApobec', sample.to_s, context.to_s, random_dataset),
          fitting_log: log_folder && File.join(log_folder, 'YeastApobec', sample.to_s, context.to_s, "#{random_dataset}.log"),
          output_file: File.join(LocalPaths::Secondary::MotifStatistics, full_folder, 'YeastApobec', sample.to_s, context.to_s, "#{random_dataset}.csv"),
          task_name: "#{task_prefix}:YeastApobec:#{sample}:#{context}:#{random_dataset}"
        )

        make_filtered_statistics_task(motif_statistics_file: File.join(LocalPaths::Secondary::MotifStatistics, full_folder, 'YeastApobec', sample.to_s, context.to_s, "#{random_dataset}.csv"),
                                      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, filtered_folder, 'YeastApobec', sample.to_s, context.to_s),
                                      task_name: "filtered_#{task_prefix}:YeastApobec:#{sample}:#{context}:#{random_dataset}")
      end

      make_all_common_motifs_tasks(
        folder: File.join(LocalPaths::Secondary::MotifStatistics, filtered_folder, 'YeastApobec', sample.to_s, context.to_s),
        output_folder: File.join(LocalPaths::Secondary::MotifStatistics, common_folder, 'YeastApobec', sample.to_s, context.to_s),
        configuration_module: Configuration::YeastApobec,
        task_name: "common_#{task_prefix}:YeastApobec:#{sample}:#{context}"
      )
    end
  end
end
