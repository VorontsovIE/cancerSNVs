desc 'Extract SNV proiles and mean information content of affected positions'
task :extract_snv_positions => ['extract_snv_positions:Alexandrov', 'extract_snv_positions:YeastApobec', 'extract_snv_positions:NikZainal']

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


def motif_infos_in_snv_positions_file(filename)
  File.readlines(filename)
      .map(&:chomp)
      .each_slice(2)
      .map{|motif_info, profile| motif_info}
end

def ic_averaged_on_core_by_motif(motif_infos)
  motif_infos.map{|motif_info|
    motif_name, ic_averaged_on_core, ic_averaged_on_everything = motif_info.split("\t").drop(1)
    [motif_name, ic_averaged_on_core.to_f]
  }.to_h
end

def ic_averaged_on_everything_by_motif(motif_infos)
  motif_infos.map{|motif_info|
    motif_name, ic_averaged_on_core, ic_averaged_on_everything = motif_info.split("\t").drop(1)
    [motif_name, ic_averaged_on_core.to_f]
  }.to_h
end

def motif_IC_matrix(cancer_ic_by_motif, random_samples_ics_by_motifs)
  random_samples = random_samples_ics_by_motifs.keys
  result = []
  result << [ 'Motif',
            'Cancer',
            *random_samples.map(&:capitalize),
            *random_samples.map{|sample| "Cancer to #{sample.downcase}"}
          ]
  motifs = cancer_ic_by_motif.keys
  motifs.each{|motif|
    cancer_ic = cancer_ic_by_motif[motif]
    random_ics = random_samples.map{|sample| random_samples_ics_by_motifs[sample][motif] }
    cancer_to_random_ics = random_ics.map{|random_ic| cancer_ic / random_ic }
    result << [ motif, cancer_ic, *random_ics, *cancer_to_random_ics ]
  }
  result
end


def samples_with_contexts_and_config
  results = []
  Configuration.getAlexandrovWholeGenomeCancers.each do |cancer_type|
    Configuration::Alexandrov.contexts_by_cancer_type(cancer_type).each do |context|
      results << [File.join('Alexandrov', cancer_type.to_s, context.to_s), Configuration::Alexandrov]
    end
  end

  Configuration.getYeastApobecSamples.each do |cancer_type|
    Configuration::YeastApobec.contexts_by_cancer_type(cancer_type).each do |context|
      results << [File.join('YeastApobec', cancer_type.to_s, context.to_s), Configuration::YeastApobec]
    end
  end

  Configuration::NikZainalContexts.each do |context|
    results << [File.join('NikZainal', context.to_s), Configuration::NikZainal]
  end
  results
end

desc 'Aggregate SNV mean KDICs'
task :aggregate_snv_positions do
  samples_with_contexts_and_config.each do |sample_path, config|
    output_folder = File.join('results/snv_positions_aggregated', sample_path)
    mkdir_p output_folder  unless Dir.exist?(output_folder)
    cancer_fn = File.join('results/snv_positions', sample_path, 'cancer.txt')
    ics_random_on_core = config.const_get(:RandomDatasets).map{|dataset|
      random_fn = File.join('results/snv_positions', sample_path, "#{dataset}.txt")
      [dataset, ic_averaged_on_core_by_motif(motif_infos_in_snv_positions_file(random_fn))]
    }.to_h
    ics_random_on_everything = config.const_get(:RandomDatasets).map{|dataset|
      random_fn = File.join('results/snv_positions', sample_path, "#{dataset}.txt")
      [dataset, ic_averaged_on_everything_by_motif(motif_infos_in_snv_positions_file(random_fn))]
    }.to_h
    motif_ics_on_core = motif_IC_matrix(
      ic_averaged_on_core_by_motif(motif_infos_in_snv_positions_file(cancer_fn)),
      ics_random_on_core
    )

    motif_ics_on_everything = motif_IC_matrix(
      ic_averaged_on_everything_by_motif(motif_infos_in_snv_positions_file(cancer_fn)),
      ics_random_on_everything
    )

    File.open(File.join(output_folder, 'mean_KDIC_on_core.tsv'), 'w') do |fw|
      print_matrix(motif_ics_on_core, stream: fw)
    end

    File.open(File.join(output_folder, 'mean_KDIC_on_everything.tsv'), 'w') do |fw|
      print_matrix(motif_ics_on_everything, stream: fw)
    end
  end
end
