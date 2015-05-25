def site_copy_task(file_from:, folder_to:, task_name: )
  directory folder_to
  to = File.join(folder_to, File.basename(file_from))
  task task_name => to
  file to => [file_from, folder_to] do
    ln file_from, to, force: true
  end
end

task 'preparations:copy_sites:Alexandrov'
AlexandrovWholeGenomeCancers.each do |cancer_type|
  Dir.glob(File.join(LocalPaths::Secondary::Chunks, 'Alexandrov', cancer_type.to_s, 'sites_*.txt')).each do |from|
    site_copy_task(file_from: from, 
                  folder_to: File.join(LocalPaths::Secondary::Sites, 'Alexandrov', cancer_type.to_s, 'any'), 
                  task_name: 'preparations:copy_sites:Alexandrov')
  end
end

task 'preparations:copy_sites:NikZainal'
Dir.glob(File.join(LocalPaths::Secondary::Chunks, 'NikZainal', 'sites_*.txt')).each do |from|
  site_copy_task(file_from: from, 
                folder_to: File.join(LocalPaths::Secondary::Sites, 'NikZainal', 'any'), 
                task_name: 'preparations:copy_sites:NikZainal')
end
