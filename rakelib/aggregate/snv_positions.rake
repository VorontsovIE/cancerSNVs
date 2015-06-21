desc 'Extract SNV proiles and mean information content of affected positions'
task :extract_snv_positions => ['extract_snv_positions:Alexandrov']

Configuration.getAlexandrovWholeGenomeCancers.each do |cancer_type|
  task 'extract_snv_positions:Alexandrov' => "extract_snv_positions:Alexandrov:#{cancer_type}"
  Configuration::Alexandrov.contexts_by_cancer_type(cancer_type).each do |context|
    task "extract_snv_positions:Alexandrov:#{cancer_type}" => "extract_snv_positions:Alexandrov:#{cancer_type}:#{context}"
    Configuration::Alexandrov::Datasets.each do |dataset|
      task "extract_snv_positions:Alexandrov:#{cancer_type}:#{context}" => "extract_snv_positions:Alexandrov:#{cancer_type}:#{context}:#{dataset}"

      output_folder = File.join('results/snv_positions', 'Alexandrov', cancer_type.to_s, context.to_s)
      directory output_folder
      task "extract_snv_positions:Alexandrov:#{cancer_type}:#{context}:#{dataset}" => [output_folder] do
        ruby 'extract_snv_position_profile.rb', 
              File.join('results/fitted_sites', 'Alexandrov', cancer_type.to_s, context.to_s, "sites_#{dataset}.txt"),
              {out: File.join(output_folder, "#{dataset}.txt")}, {}
      end
    end
  end
end

Configuration.getYeastApobecSamples.each do |cancer_type|
  task 'extract_snv_positions:YeastApobec' => "extract_snv_positions:YeastApobec:#{cancer_type}"
  Configuration::YeastApobec.contexts_by_cancer_type(cancer_type).each do |context|
    task "extract_snv_positions:YeastApobec:#{cancer_type}" => "extract_snv_positions:YeastApobec:#{cancer_type}:#{context}"
    Configuration::YeastApobec::Datasets.each do |dataset|
      task "extract_snv_positions:YeastApobec:#{cancer_type}:#{context}" => "extract_snv_positions:YeastApobec:#{cancer_type}:#{context}:#{dataset}"

      output_folder = File.join('results/snv_positions', 'YeastApobec', cancer_type.to_s, context.to_s)
      directory output_folder
      task "extract_snv_positions:YeastApobec:#{cancer_type}:#{context}:#{dataset}" => [output_folder] do
        ruby 'extract_snv_position_profile.rb', 
              File.join('results/fitted_sites', 'YeastApobec', cancer_type.to_s, context.to_s, "sites_#{dataset}.txt"),
              {out: File.join(output_folder, "#{dataset}.txt")}, {}
      end
    end
  end
end

Configuration::NikZainalContexts.each do |context|
  task "extract_snv_positions:NikZainal" => "extract_snv_positions:NikZainal:#{context}"
  Configuration::NikZainal::Datasets.each do |dataset|
    task "extract_snv_positions:NikZainal:#{context}" => "extract_snv_positions:NikZainal:#{context}:#{dataset}"

    output_folder = File.join('results/snv_positions', 'NikZainal', context.to_s)
    directory output_folder
    task "extract_snv_positions:NikZainal:#{context}:#{dataset}" => [output_folder] do
      ruby 'extract_snv_position_profile.rb', 
            File.join('results/fitted_sites', 'NikZainal', context.to_s, "sites_#{dataset}.txt"),
            {out: File.join(output_folder, "#{dataset}.txt")}, {}
    end
  end
end
