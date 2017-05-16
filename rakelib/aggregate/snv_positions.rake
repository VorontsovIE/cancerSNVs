desc 'Extract SNV profiles and mean information content of affected positions'
task :extract_snv_positions

Configuration.WholeGenomeCancers.each do |cancer_type|
  task 'extract_snv_positions' => "extract_snv_positions:#{cancer_type}"
  Configuration::Datasets.each do |dataset|
    task "extract_snv_positions:#{cancer_type}" => "extract_snv_positions:#{cancer_type}:#{dataset}"

    output_folder = File.join(LocalPaths::Results, 'snv_positions', cancer_type.to_s)
    directory output_folder
    task "extract_snv_positions:#{cancer_type}:#{dataset}" => [output_folder] do
      ruby 'extract_snv_position_profile.rb',
            File.join(LocalPaths::Results, 'fitted_sites', cancer_type.to_s, "sites_#{dataset}.txt"),
            {out: File.join(output_folder, "#{dataset}.txt")}, {}
    end
  end
end
