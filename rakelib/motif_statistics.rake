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
                        '--expand-control-set', Configuration::ExpandControlSetFold,
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
        ruby 'common_motifs.rb', *input_files, {out: output_file}, {}
      end
    end
  end
end

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
    task "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs" => [
      "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs:genome",
      "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs:shuffle",
      "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs:all"
    ]


    make_common_motifs_tasks(
      folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'Alexandrov', cancer_type.to_s, context.to_s),
      random_datasets: Configuration::Alexandrov::RandomGenomeDatasets,
      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'Alexandrov', cancer_type.to_s, context.to_s),
      output_file: 'compared_to_each_genome.txt',
      task_name: "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs:genome"
    )

    make_common_motifs_tasks(
      folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'Alexandrov', cancer_type.to_s, context.to_s),
      random_datasets: Configuration::Alexandrov::RandomShuffleDatasets,
      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'Alexandrov', cancer_type.to_s, context.to_s),
      output_file: 'compared_to_each_shuffle.txt',
      task_name: "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs:shuffle"
    )


    make_common_motifs_tasks(
      folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'Alexandrov', cancer_type.to_s, context.to_s),
      random_datasets: Configuration::Alexandrov::RandomDatasets,
      output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'Alexandrov', cancer_type.to_s, context.to_s),
      output_file: 'compared_to_each.txt',
      task_name: "motif_statistics:Alexandrov:#{cancer_type}:#{context}:common_motifs:all"
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
  task "motif_statistics:NikZainal:#{context}:common_motifs" => [
    "motif_statistics:NikZainal:#{context}:common_motifs:genome",
    "motif_statistics:NikZainal:#{context}:common_motifs:shuffle",
    "motif_statistics:NikZainal:#{context}:common_motifs:all"
  ]


  make_common_motifs_tasks(
    folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'NikZainal', context.to_s),
    random_datasets: Configuration::NikZainal::RandomGenomeDatasets,
    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'NikZainal', context.to_s),
    output_file: 'compared_to_each_genome.txt',
    task_name: "motif_statistics:NikZainal:#{context}:common_motifs:genome"
  )

  make_common_motifs_tasks(
    folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'NikZainal', context.to_s),
    random_datasets: Configuration::NikZainal::RandomShuffleDatasets,
    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'NikZainal', context.to_s),
    output_file: 'compared_to_each_shuffle.txt',
    task_name: "motif_statistics:NikZainal:#{context}:common_motifs:shuffle"
  )


  make_common_motifs_tasks(
    folder: File.join(LocalPaths::Secondary::MotifStatistics, 'filtered', 'NikZainal', context.to_s),
    random_datasets: Configuration::NikZainal::RandomDatasets,
    output_folder: File.join(LocalPaths::Secondary::MotifStatistics, 'common', 'NikZainal', context.to_s),
    output_file: 'compared_to_each.txt',
    task_name: "motif_statistics:NikZainal:#{context}:common_motifs:all"
  )
end
