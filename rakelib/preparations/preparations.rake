namespace :preparations do
  desc 'Prepare all data (SNVs: cancer and random and generate chunks to run on multiple cores)'
  task :all do
    Rake::Task['preparations:extractSNVs'].invoke
    Rake::Task['preparations:generate_random_SNVs'].invoke
    Rake::Task['preparations:generate_chunks'].invoke
    Rake::Task['preparations:shuffle_sites'].invoke
  end
end
