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
task 'slicing'

WholeGenomeCancers.each do |cancer_type|
  task 'slicing' => "slicing:#{cancer_type}"
  
  Configuration::Datasets.each do |dataset|
    task "slicing:#{cancer_type}" => "slicing:#{cancer_type}:#{dataset}"

    make_slicing_task(sites_file: File.join(LocalPaths::Secondary::Fitting, cancer_type.to_s, "sites_#{dataset}.txt"),
                      output_folder: File.join(LocalPaths::Secondary::Slices, cancer_type.to_s, dataset),
                      task_name: "slicing:#{cancer_type}:#{dataset}")
  end
end
