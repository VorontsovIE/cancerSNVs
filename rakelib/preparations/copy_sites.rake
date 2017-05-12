def site_shuffle_task(file_from:, folder_to:, task_name: )
  directory folder_to
  file_to = File.join(folder_to, File.basename(file_from))
  task task_name => file_to
  file file_to => [file_from, folder_to] do
    sh 'shuf', '--output', file_to, file_from
  end
end

desc 'Shuffle sites from perfectosape output (from chunks folder) and put them to sites folder'
task 'preparations:shuffle_sites' => ['preparations:shuffle_sites:Alexandrov']

task 'preparations:shuffle_sites:Alexandrov'
AlexandrovWholeGenomeCancers.each do |cancer_type|
  Dir.glob(File.join(LocalPaths::Secondary::Chunks, 'Alexandrov', cancer_type.to_s, 'sites_*.txt')).each do |from|
    site_shuffle_task(file_from: from,
                  folder_to: File.join(LocalPaths::Secondary::Sites, 'Alexandrov', cancer_type.to_s, 'any'),
                  task_name: 'preparations:shuffle_sites:Alexandrov')
  end
end
