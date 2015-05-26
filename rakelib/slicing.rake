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
task 'slicing' => ['slicing:Alexandrov', 'slicing:NikZainal']

task 'slicing:Alexandrov'
AlexandrovWholeGenomeCancers.each do |cancer_type|
  task 'slicing:Alexandrov' => "slicing:Alexandrov:#{cancer_type}"
  Configuration::Alexandrov.contexts_by_cancer_type(cancer_type).each do |context|
    task "slicing:Alexandrov:#{cancer_type}" => "slicing:Alexandrov:#{cancer_type}:#{context}"

    Configuration::Alexandrov::Datasets.each do |dataset|
      task "slicing:Alexandrov:#{cancer_type}:#{context}" => "slicing:Alexandrov:#{cancer_type}:#{context}:#{dataset}"

      make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, 'Alexandrov', cancer_type.to_s, context.to_s, "sites_#{dataset}.txt"),
                        output_folder: File.join(LocalPaths::Secondary::Slices, 'Alexandrov', cancer_type.to_s, context.to_s, dataset),
                        task_name: "slicing:Alexandrov:#{cancer_type}:#{context}:#{dataset}")
    end
  end
end

task 'slicing:NikZainal'
Configuration::NikZainalContexts.each do |context|
  task 'slicing:NikZainal' => "slicing:NikZainal:#{context}"
  Configuration::NikZainal::Datasets.each do |dataset|
    task "slicing:NikZainal:#{context}" => "slicing:NikZainal:#{context}:#{dataset}"

    make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, 'NikZainal', context.to_s, "sites_#{dataset}.txt"),
                      output_folder: File.join(LocalPaths::Secondary::Slices, 'NikZainal', context.to_s, dataset),
                      task_name: "slicing:NikZainal:#{context}:#{dataset}")
  end
end
