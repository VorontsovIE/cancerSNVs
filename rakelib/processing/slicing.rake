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
task 'slicing' => ['slicing:Alexandrov']
desc 'Make slices of motif statistics (w/o fitting)'
task 'slicing_wo_fitting' => ['slicing_wo_fitting:Alexandrov']

desc 'Make slices of motif statistics for Alexandrov datasets'
task 'slicing:Alexandrov'
desc 'Make slices of motif statistics for Alexandrov datasets (w/o fitting)'
task 'slicing_wo_fitting:Alexandrov'

AlexandrovWholeGenomeCancers.each do |cancer_type|
  task 'slicing:Alexandrov' => "slicing:Alexandrov:#{cancer_type}"
  task 'slicing_wo_fitting:Alexandrov' => "slicing_wo_fitting:Alexandrov:#{cancer_type}"
  
  Configuration::Alexandrov::Datasets.each do |dataset|
    task "slicing:Alexandrov:#{cancer_type}" => "slicing:Alexandrov:#{cancer_type}:#{dataset}"
    task "slicing_wo_fitting:Alexandrov:#{cancer_type}" => "slicing_wo_fitting:Alexandrov:#{cancer_type}:#{dataset}"

    make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, 'Alexandrov', cancer_type.to_s, "sites_#{dataset}.txt"),
                      output_folder: File.join(LocalPaths::Secondary::Slices, 'Alexandrov', cancer_type.to_s, dataset),
                      task_name: "slicing:Alexandrov:#{cancer_type}:#{dataset}")

    make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Sites, 'Alexandrov', cancer_type.to_s, "sites_#{dataset}.txt"),
                      output_folder: File.join(LocalPaths::Results, 'motif_statistics/slices_wo_fitting', 'Alexandrov', cancer_type.to_s, dataset),
                      task_name: "slicing_wo_fitting:Alexandrov:#{cancer_type}:#{dataset}")
  end
end
