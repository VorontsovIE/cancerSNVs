def make_slicing_task(sites_file:, output_folder:, task_name:)
  directory output_folder
  task task_name => [output_folder, sites_file] do
    ruby 'calculate_motif_statistics.rb',
        sites_file, output_folder,
        '--pvalue-cutoff', Configuration::PvalueCutoff.to_s,
        '--fold-change-cutoff', Configuration::FoldChangeCutoff.to_s
  end
end

desc 'Make slices of motif statistics'
task 'slicing' => ['slicing:Alexandrov', 'slicing:NikZainal', 'slicing:YeastApobec']
desc 'Make slices of motif statistics (w/o fitting)'
task 'slicing_wo_fitting' => ['slicing_wo_fitting:Alexandrov', 'slicing_wo_fitting:NikZainal', 'slicing_wo_fitting:YeastApobec']

desc 'Make slices of motif statistics for Alexandrov datasets'
task 'slicing:Alexandrov'
desc 'Make slices of motif statistics for Alexandrov datasets (w/o fitting)'
task 'slicing_wo_fitting:Alexandrov'

AlexandrovWholeGenomeCancers.each do |cancer_type|
  task 'slicing:Alexandrov' => "slicing:Alexandrov:#{cancer_type}"
  task 'slicing_wo_fitting:Alexandrov' => "slicing_wo_fitting:Alexandrov:#{cancer_type}"
  Configuration::Alexandrov.contexts_by_cancer_type(cancer_type).each do |context|
    task "slicing:Alexandrov:#{cancer_type}" => "slicing:Alexandrov:#{cancer_type}:#{context}"
    task "slicing_wo_fitting:Alexandrov:#{cancer_type}" => "slicing_wo_fitting:Alexandrov:#{cancer_type}:#{context}"

    Configuration::Alexandrov::Datasets.each do |dataset|
      task "slicing:Alexandrov:#{cancer_type}:#{context}" => "slicing:Alexandrov:#{cancer_type}:#{context}:#{dataset}"
      task "slicing_wo_fitting:Alexandrov:#{cancer_type}:#{context}" => "slicing_wo_fitting:Alexandrov:#{cancer_type}:#{context}:#{dataset}"

      make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, 'Alexandrov', cancer_type.to_s, context.to_s, "sites_#{dataset}.txt"),
                        output_folder: File.join(LocalPaths::Secondary::Slices, 'Alexandrov', cancer_type.to_s, context.to_s, dataset),
                        task_name: "slicing:Alexandrov:#{cancer_type}:#{context}:#{dataset}")

      make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Sites, 'Alexandrov', cancer_type.to_s, context.to_s, "sites_#{dataset}.txt"),
                        output_folder: File.join(LocalPaths::Results, 'motif_statistics/slices_wo_fitting', 'Alexandrov', cancer_type.to_s, context.to_s, dataset),
                        task_name: "slicing_wo_fitting:Alexandrov:#{cancer_type}:#{context}:#{dataset}")
    end
  end
end

desc 'Make slices of motif statistics for Nik-Zainal dataset'
task 'slicing:NikZainal'
desc 'Make slices of motif statistics for Nik-Zainal dataset (w/o fitting)'
task 'slicing_wo_fitting:NikZainal'

Configuration::NikZainalContexts.each do |context|
  task 'slicing:NikZainal' => "slicing:NikZainal:#{context}"
  task 'slicing_wo_fitting:NikZainal' => "slicing_wo_fitting:NikZainal:#{context}"
  Configuration::NikZainal::Datasets.each do |dataset|
    task "slicing:NikZainal:#{context}" => "slicing:NikZainal:#{context}:#{dataset}"
    task "slicing_wo_fitting:NikZainal:#{context}" => "slicing_wo_fitting:NikZainal:#{context}:#{dataset}"

    make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, 'NikZainal', context.to_s, "sites_#{dataset}.txt"),
                      output_folder: File.join(LocalPaths::Secondary::Slices, 'NikZainal', context.to_s, dataset),
                      task_name: "slicing:NikZainal:#{context}:#{dataset}")

    make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Sites, 'NikZainal', context.to_s, "sites_#{dataset}.txt"),
                      output_folder: File.join(LocalPaths::Results, 'motif_statistics/slices_wo_fitting', 'NikZainal', context.to_s, dataset),
                      task_name: "slicing_wo_fitting:NikZainal:#{context}:#{dataset}")
  end
end

desc 'Make slices of motif statistics for APOBEC mutations in yeast'
task 'slicing:YeastApobec'
desc 'Make slices of motif statistics for APOBEC mutations in yeast (w/o fitting)'
task 'slicing_wo_fitting:YeastApobec'

YeastApobecSamples.each do |sample|
  task 'slicing:YeastApobec' => "slicing:YeastApobec:#{sample}"
  task 'slicing_wo_fitting:YeastApobec' => "slicing_wo_fitting:YeastApobec:#{sample}"
  Configuration::YeastApobec.contexts_by_cancer_type(sample).each do |context|
    task "slicing:YeastApobec:#{sample}" => "slicing:YeastApobec:#{sample}:#{context}"
    task "slicing_wo_fitting:YeastApobec:#{sample}" => "slicing_wo_fitting:YeastApobec:#{sample}:#{context}"

    Configuration::YeastApobec::Datasets.each do |dataset|
      task "slicing:YeastApobec:#{sample}:#{context}" => "slicing:YeastApobec:#{sample}:#{context}:#{dataset}"
      task "slicing_wo_fitting:YeastApobec:#{sample}:#{context}" => "slicing_wo_fitting:YeastApobec:#{sample}:#{context}:#{dataset}"

      make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, 'YeastApobec', sample.to_s, context.to_s, "sites_#{dataset}.txt"),
                        output_folder: File.join(LocalPaths::Secondary::Slices, 'YeastApobec', sample.to_s, context.to_s, dataset),
                        task_name: "slicing:YeastApobec:#{sample}:#{context}:#{dataset}")

      make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Sites, 'YeastApobec', sample.to_s, context.to_s, "sites_#{dataset}.txt"),
                        output_folder: File.join(LocalPaths::Results, 'motif_statistics/slices_wo_fitting', 'YeastApobec', sample.to_s, context.to_s, dataset),
                        task_name: "slicing_wo_fitting:YeastApobec:#{sample}:#{context}:#{dataset}")
    end
  end
end
