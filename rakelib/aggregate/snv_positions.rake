require 'bioinform'
require_relative '../../lib/perfectosape/results_short'

def extract_snv_positions_profile(sites_filename)
  results = Hash.new{|hsh, motif| hsh[motif] = Hash.new(0) }
  PerfectosAPE::Result.each_in_file(sites_filename) do |site_info|
    pos = (site_info.pvalue_1 <= site_info.pvalue_2) ? site_info.snv_position_in_site_1_pwm :  site_info.snv_position_in_site_2_pwm
    results[site_info.motif_name][pos] += 1
  end
  results
end

def dic(pos)
  n = pos.inject(0.0, &:+)
  (pos.map{|el| Math.lgamma(1 + el).first }.inject(0.0, &:+) - Math.lgamma(1 + n).first) * 1.0 / n
end

def kdic(pos, background: [0.25, 0.25, 0.25, 0.25])
  n = pos.inject(0.0, &:+)
  dic(pos) - pos.each_with_index.map{|el, ind| el * Math.log(background[ind]) / n }.inject(0.0, &:+)
end


def extract_snv_positions_profile_to_file(sites_filename, output_filename, motif_pcms, expansion_flank_length: 11)
  profiles = extract_snv_positions_profile(sites_filename)
  File.open(output_filename, 'w') do |fw|
    profiles.sort.each{|motif, profile_hash|
      motif_pcm = motif_pcms[motif]
      total_substitutions = profile_hash.each_value.inject(0, &:+)
      total_substitutions_in_core = (0 ... motif_pcm.length).map{|pos| profile_hash[pos] }.inject(0, &:+)

      total_dic = (0 ... motif_pcm.length).map{|pos|
        kdic(motif_pcm.matrix[pos]) * profile_hash[pos]
      }.inject(0.0, &:+)

      profile = (-expansion_flank_length ... (motif_pcms[motif].length + expansion_flank_length)).map{|pos|
        profile_hash[pos]
      }

      fw.puts "> #{motif}"
      fw.puts profile.join("\t")
      fw.puts (total_dic.to_f / total_substitutions_in_core)
      fw.puts (total_dic.to_f / total_substitutions)
    }
  end
end

motif_pcms = {}
Dir.glob(File.join(LocalPaths::MotifCollectionPCM, '*.pcm')) do |motif_filename|
  motif = Bioinform::MotifModel::PCM.from_file(motif_filename)
  motif_pcms[motif.name.to_sym] ||= motif
end

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
        sites_filename = File.join('results/fitted_sites', 'Alexandrov', cancer_type.to_s, context.to_s, "#{dataset}.txt")
        output_filename = File.join(output_folder, "#{dataset}.txt")
        extract_snv_positions_profile_to_file(sites_filename, output_filename, motif_pcms, expansion_flank_length: Configuration::ExpandFlanksLength)
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
        sites_filename = File.join('results/fitted_sites', 'YeastApobec', cancer_type.to_s, context.to_s, "#{dataset}.txt")
        output_filename = File.join(output_folder, "#{dataset}.txt")
        extract_snv_positions_profile_to_file(sites_filename, output_filename, motif_pcms, expansion_flank_length: Configuration::ExpandFlanksLength)
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
      sites_filename = File.join('results/fitted_sites', 'NikZainal', context.to_s, "#{dataset}.txt")
      output_filename = File.join(output_folder, "#{dataset}.txt")
      extract_snv_positions_profile_to_file(sites_filename, output_filename, motif_pcms, expansion_flank_length: Configuration::ExpandFlanksLength)
    end
  end
end
