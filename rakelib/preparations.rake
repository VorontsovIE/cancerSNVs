namespace :preparations do
  task :all do
    Rake::Task['preparations:extractSNVs'].invoke
    Rake::Task['preparations:generate_random_SNVs'].invoke
    Rake::Task['preparations:generate_chunks'].invoke
  end
end
