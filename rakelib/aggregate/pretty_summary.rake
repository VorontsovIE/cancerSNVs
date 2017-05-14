desc 'Make summary for all experiments in a single table'
task 'pretty_summary' do
  SUMMARY_HEADERS = {
    motif: 'Motif name (HOCOMOCO v9)',
    gene: 'TF gene (Gene symbol)',
    tf_family: 'TF family (TFClass)',

    affinity_loss_sites: 'Predicted TFBS, affinity loss cases',
    germline_allele_sites: 'Predicted TFBS, total for germline alleles',
    affinity_loss_sites_genomic_control: 'Control (genomic) predicted TFBS, affinity loss cases',
    germline_allele_sites_genomic_control: 'Control (genomic) predicted TFBS, total for germline alleles',
    affinity_loss_sites_shuffle_control: 'Control (shuffle) predicted TFBS, affinity loss cases',
    germline_allele_sites_shuffle_control: 'Control (shuffle) predicted TFBS, total for germline alleles',

    affinity_gain_sites: 'Predicted TFBS, affinity gain cases',
    mutated_allele_sites: 'Predicted TFBS, total for mutated alleles',
    affinity_gain_sites_genomic_control: 'Control (genomic) predicted TFBS, affinity gain cases',
    mutated_allele_sites_genomic_control: 'Control (genomic) predicted TFBS, total for induced substitutions',
    affinity_gain_sites_shuffle_control: 'Control (shuffle) predicted TFBS, affinity gain cases',
    mutated_allele_sites_shuffle_control: 'Control (shuffle) predicted TFBS, total for induced substitutions',

    mutations_within_site: 'Mutations within TFBS, germline alleles',
    mutations_within_site_genomic_control: 'Control (genomic), substitutions within TFBS, germline alleles',
    mutations_within_site_shuffle_control: 'Control (shuffle), substitutions within TFBS, germline alleles',
    

    affinity_loss_rate_genomic_control: 'Affinity loss, relative rate (cancer/genomic control)',
    affinity_loss_rate_shuffle_control: 'Affinity loss, relative rate (cancer/shuffle control)',

    affinity_gain_rate_genomic_control: 'Affinity gain, relative rate (cancer/genomic control)',
    affinity_gain_rate_shuffle_control: 'Affinity gain, relative rate (cancer/shuffle control)',
    
    mutations_within_site_relative_rate_genomic_control: 'Mutations within motif, relative frequency (cancer/genomic control)',
    mutations_within_site_relative_rate_shuffle_control: 'Mutations within motif, relative frequency (cancer/shuffle control)',


    affinity_loss_significance_genomic_control: "Affinity loss significance (FDR-corrected two-tail Fisher's exact test P-value), versus genomic control",
    affinity_loss_significance_shuffle_control: "Affinity loss significance (FDR-corrected two-tail Fisher's exact test P-value), versus shuffle control",

    affinity_gain_significance_genomic_control: "Affinity gain significance (FDR-corrected two-tail Fisher's exact test P-value), versus genomic control",
    affinity_gain_significance_shuffle_control: "Affinity gain significance (FDR-corrected two-tail Fisher's exact test P-value), versus shuffle control",

    mutations_within_site_significance_genomic_control: "Significance (FDR-corrected two-tail Fisher's exact test P-value), versus genomic control",
    mutations_within_site_significance_shuffle_control: "Significance (FDR-corrected two-tail Fisher's exact test P-value), versus shuffle control",

    underfitted_sites_genomic_control: 'TFBS missing from genomic control',
    underfitted_sites_shuffle_control: 'TFBS missing from shuffle control',
  }

  affinity_infos_header = [
    :motif,
    :gene,
    :tf_family,

    :affinity_loss_sites,
    :germline_allele_sites,
    :affinity_loss_sites_genomic_control,
    :germline_allele_sites_genomic_control,
    :affinity_loss_sites_shuffle_control,
    :germline_allele_sites_shuffle_control,

    :affinity_gain_sites,
    :mutated_allele_sites,
    :affinity_gain_sites_genomic_control,
    :mutated_allele_sites_genomic_control,
    :affinity_gain_sites_shuffle_control,
    :mutated_allele_sites_shuffle_control,

    :affinity_loss_rate_genomic_control,
    :affinity_loss_rate_shuffle_control,

    :affinity_gain_rate_genomic_control,
    :affinity_gain_rate_shuffle_control,
    
    :affinity_loss_significance_genomic_control,
    :affinity_loss_significance_shuffle_control,

    :affinity_gain_significance_genomic_control,
    :affinity_gain_significance_shuffle_control,
    :underfitted_sites_genomic_control,
    :underfitted_sites_shuffle_control,
  ].map{|column| SUMMARY_HEADERS[column] }

  position_infos_header = [
    :motif,
    :gene,
    :tf_family,

    :mutations_within_site,
    :germline_allele_sites,
    :mutations_within_site_genomic_control,
    :germline_allele_sites_genomic_control,
    :mutations_within_site_shuffle_control,
    :germline_allele_sites_shuffle_control,
    
    :mutations_within_site_relative_rate_genomic_control,
    :mutations_within_site_relative_rate_shuffle_control,

    :mutations_within_site_significance_genomic_control,
    :mutations_within_site_significance_shuffle_control,
    :underfitted_sites_genomic_control,
    :underfitted_sites_shuffle_control,
  ].map{|column| SUMMARY_HEADERS[column] }

  Configuration.WholeGenomeCancers.each do |cancer_type|
    folder = File.join(LocalPaths::Results, "motif_statistics/full/#{cancer_type}/")
    output_folder = File.join(LocalPaths::Results, "motif_statistics/pretty_summary/")

    genome_control_infos = File.readlines(File.join(folder, "random_genome.csv")).drop(1)
    shuffle_control_infos = File.readlines(File.join(folder, "random_shuffle.csv")).drop(1)

    infos = genome_control_infos.zip(shuffle_control_infos).map{|genome_control_line, shuffle_control_line|
      motif, gene, quality, \
        affinity_loss_rate_genomic_control, affinity_loss_significance_genomic_control, \
        affinity_gain_rate_genomic_control, affinity_gain_significance_genomic_control, \
        mutations_within_site_relative_rate_genomic_control, mutations_within_site_significance_genomic_control, \

        affinity_loss_sites_genomic_control, \
        affinity_gain_sites_genomic_control, \
        germline_allele_sites_genomic_control, \
        mutated_allele_sites_genomic_control, \
        mutations_within_site_genomic_control, \
        mutations_outside_site_genomic_control, \
      
        affinity_loss_sites, \
        affinity_gain_sites, \
        germline_allele_sites, \
        mutated_allele_sites, \
        mutations_within_site, \
        mutations_outside_site, \
        *rest =  genome_control_line.split("\t")

        underfitted_sites_genomic_control = rest[10]

      #########
      _motif, _gene, _quality, \
        affinity_loss_rate_shuffle_control, affinity_loss_significance_shuffle_control, \
        affinity_gain_rate_shuffle_control, affinity_gain_significance_shuffle_control, \
        mutations_within_site_relative_rate_shuffle_control, mutations_within_site_significance_shuffle_control, \

        affinity_loss_sites_shuffle_control, \
        affinity_gain_sites_shuffle_control, \
        germline_allele_sites_shuffle_control, \
        mutated_allele_sites_shuffle_control, \
        mutations_within_site_shuffle_control, \
        mutations_outside_site_shuffle_control, \
      
        _affinity_loss_sites, \
        _affinity_gain_sites, \
        _germline_allele_sites, \
        _mutated_allele_sites, \
        _mutations_within_site, \
        _mutations_outside_site, \
        *rest =  shuffle_control_line.split("\t")

        underfitted_sites_shuffle_control = rest[10]

      #########
      common_data_1 = [
        motif, gene, quality,
        affinity_loss_sites, affinity_gain_sites, 
        germline_allele_sites, mutated_allele_sites, 
        mutations_within_site, mutations_outside_site
      ]
      common_data_2 = [
        _motif, _gene, _quality,
        _affinity_loss_sites, _affinity_gain_sites, 
        _germline_allele_sites, _mutated_allele_sites, 
        _mutations_within_site, _mutations_outside_site
      ]

      raise 'Data in tables which should be the same differs'  unless common_data_1 == common_data_2
      tf_family = MOTIF_FAMILY_RECOGNIZERS[3].subfamilies_by_motif(motif).map(&:to_s).join('; ')

      # # Don't recalculate them!
      # mutations_within_site_rate = (mutations_within_site.to_f / germline_allele_sites.to_f)
      # mutations_within_site_rate_genomic_control = (mutations_within_site_genomic_control.to_f / germline_allele_sites_genomic_control.to_f)
      # mutations_within_site_rate_shuffle_control = (mutations_within_site_shuffle_control.to_f / germline_allele_sites_shuffle_control.to_f)

      # mutations_within_site_relative_rate_genomic_control = mutations_within_site_rate_genomic_control / mutations_within_site_rate
      # mutations_within_site_relative_rate_shuffle_control = mutations_within_site_rate_shuffle_control / mutations_within_site_rate

      affinity_infos = [
        motif,
        gene,
        tf_family,

        affinity_loss_sites,
        germline_allele_sites,
        affinity_loss_sites_genomic_control,
        germline_allele_sites_genomic_control,
        affinity_loss_sites_shuffle_control,
        germline_allele_sites_shuffle_control,

        affinity_gain_sites,
        mutated_allele_sites,
        affinity_gain_sites_genomic_control,
        mutated_allele_sites_genomic_control,
        affinity_gain_sites_shuffle_control,
        mutated_allele_sites_shuffle_control,

        affinity_loss_rate_genomic_control,
        affinity_loss_rate_shuffle_control,

        affinity_gain_rate_genomic_control,
        affinity_gain_rate_shuffle_control,
        
        affinity_loss_significance_genomic_control,
        affinity_loss_significance_shuffle_control,

        affinity_gain_significance_genomic_control,
        affinity_gain_significance_shuffle_control,

        underfitted_sites_genomic_control,
        underfitted_sites_shuffle_control,
      ]

      position_infos = [
        motif,
        gene,
        tf_family,

        mutations_within_site,
        germline_allele_sites,
        mutations_within_site_genomic_control,
        germline_allele_sites_genomic_control,
        mutations_within_site_shuffle_control,
        germline_allele_sites_shuffle_control,
        
        mutations_within_site_relative_rate_genomic_control,
        mutations_within_site_relative_rate_shuffle_control,

        mutations_within_site_significance_genomic_control,
        mutations_within_site_significance_shuffle_control,

        underfitted_sites_genomic_control,
        underfitted_sites_shuffle_control,
      ]

      [affinity_infos, position_infos]
    }

    FileUtils.mkdir_p File.join(output_folder, 'affinity_gain_loss_test')
    FileUtils.mkdir_p File.join(output_folder, 'core_flank_test')

    File.open(File.join(output_folder, 'affinity_gain_loss_test', "#{cancer_type}.tsv"), 'w') do |fw|
      fw.puts affinity_infos_header.join("\t")
      infos.each do |affinity_infos, position_infos|
        fw.puts affinity_infos.join("\t")
      end
    end
    
    File.open(File.join(output_folder, 'core_flank_test', "#{cancer_type}.tsv"), 'w') do |fw|
      fw.puts position_infos_header.join("\t")
      infos.each do |affinity_infos, position_infos|
        fw.puts position_infos.join("\t")
      end
    end
  end
end
