require 'bioinform'
require_relative 'experiment_configuration'
require_relative 'lib/perfectosape/results_short'

def extract_snv_positions_profile(sites_filename)
  results = Hash.new{|hsh, motif| hsh[motif] = Hash.new(0) }
  PerfectosAPE::ResultShort.each_in_file(sites_filename) do |site_info|
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


def print_snv_positions_profile(sites_filename, motif_pcms, expansion_flank_length: 11)
  profiles = extract_snv_positions_profile(sites_filename)

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

    puts ['>', motif, (total_dic.to_f / total_substitutions_in_core), (total_dic.to_f / total_substitutions)].join("\t")
    puts profile.join("\t")
  }
end

motif_pcms = {}
Dir.glob(File.join(LocalPaths::MotifCollectionPCM, '*.pcm')) do |motif_filename|
  motif = Bioinform::MotifModel::PCM.from_file(motif_filename)
  motif_pcms[motif.name.to_sym] ||= motif
end

sites_filename = ARGV[0]
print_snv_positions_profile(sites_filename, motif_pcms, expansion_flank_length: Configuration::ExpandFlanksLength)
