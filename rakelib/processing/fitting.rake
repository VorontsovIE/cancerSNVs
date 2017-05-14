def copy_cancer_sites_as_fitted_task(folder_from:, folder_to:, task_name:)
  directory folder_to
  task task_name => File.join(folder_to, 'sites_cancer.txt')
  file File.join(folder_to, 'sites_cancer.txt') => [File.join(folder_from, 'sites_cancer.txt'), folder_to] do
    ln File.join(folder_from, 'sites_cancer.txt'), File.join(folder_to, 'sites_cancer.txt'), force: true
  end
end

# fit random sites
def fit_sites_task(random_dataset:, fold:, folder_from:, folder_to:, log_folder:, task_name:)
  sites_cancer_fn = File.join(folder_from, 'sites_cancer.txt')
  sites_random_fn = File.join(folder_from, "sites_#{random_dataset}.txt")

  output_file = File.join(folder_to, "sites_#{random_dataset}.txt")
  log_file = File.join(log_folder, "#{random_dataset}.log")

  directory folder_to
  directory log_folder

  task task_name => output_file
  # should we specify log-file as output too or it'll cause duplication of actions?
  file output_file => [sites_cancer_fn, sites_random_fn, folder_to, log_folder] do
    ruby  'fitting_random_sites.rb',
          sites_cancer_fn, sites_random_fn,
          '--fold', fold.to_s,
          '--pvalue-cutoff', Configuration::PvalueCutoff.to_s,
          { out: output_file,
            err: log_file },
          {}
  end
end

desc 'Fit random sites to make control datasets with same site and mutation context rates as cancer ones'
task 'fitting' => ['fitting:Alexandrov']

task 'fitting:Alexandrov'
WholeGenomeCancers.each do |cancer_type|
  task 'fitting:Alexandrov' => "fitting:Alexandrov:#{cancer_type}"
  
  folder_from = File.join(LocalPaths::Secondary::Sites, 'Alexandrov', cancer_type.to_s)
  folder_to = File.join(LocalPaths::Secondary::Fitting, 'Alexandrov', cancer_type.to_s)

  copy_cancer_sites_as_fitted_task(folder_from: folder_from, folder_to: folder_to, task_name: "fitting:Alexandrov:#{cancer_type}")

  Configuration::Alexandrov::RandomDatasets.each do |random_dataset|
    case random_dataset
    when /genome/
      fold = Configuration::Alexandrov::FittingFoldGenome[cancer_type]
    when /shuffle/
      fold = Configuration::Alexandrov::FittingFoldShuffle[cancer_type]
    else
      next
    end
    fit_sites_task(random_dataset: random_dataset, fold: fold,
                  folder_from: folder_from, folder_to: folder_to,
                  log_folder: File.join(LocalPaths::Secondary::LogFolder, 'Alexandrov', cancer_type.to_s),
                  task_name: "fitting:Alexandrov:#{cancer_type}")
  end
end
