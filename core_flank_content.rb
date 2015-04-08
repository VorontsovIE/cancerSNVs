$:.unshift File.absolute_path('./lib', __dir__)
require 'site_info'
require 'bioinform'
require 'fisher_table'
require 'statistical_significance'
require 'optparse'
# require 'rate_comparison_infos'

motif_collection_folder = './source_data/motif_collection/'

motifs = Dir.glob(File.join(motif_collection_folder, '*.pwm')).map{|fn| Bioinform::MotifModel::PWM.from_file(fn) }
motif_names = motifs.map(&:name).map(&:to_sym)

class MutatatedSiteInfo
  def mutation_in_core?
    # raise 'Possibly you use a procedure appliable for extended motifs on unextended ones'  if length % 2 != 0
    # raise 'Possibly you use a procedure appliable for extended motifs on unextended ones'  if length <= 22
    pos = snv_position_in_site_1_pwm
    # mot_core_length = length / 2
    # (pos >= mot_core_length / 2) && pos < mot_core_length + mot_core_length / 2
    pos >= 0 && pos < length
  end
  def mutation_in_flank?
    ! in_core?
  end
end

MotifStat = Struct.new(:cancer_in_core, :random_in_core, :cancer_in_flank, :random_in_flank) do
  def self.empty
    self.new(0, 0, 0, 0)
  end

  def cancer_core_rate
    fisher_table.class_a_positive_rate
  end

  def random_core_rate
    fisher_table.class_b_positive_rate
  end

  def cancer_to_random_rate
    fisher_table.a_to_b_positive_rate_ratio
  end

  def fisher_table
    FisherTable.by_two_classes(class_a_positive: cancer_in_core, class_b_positive: random_in_core,
                                class_a_negative: cancer_in_flank, class_b_negative: random_in_flank)
  end

  def significance
    fisher_table.significance
  end
end

significance_rate_cutoff = 1.0
OptionParser.new{|opts|
  opts.on('--significance CUTOFF', 'Maximal significance rate') {|value|
    significance_rate_cutoff = Float(value)
  }
}.parse!(ARGV)

raise 'Specify cancer site infos file'  unless cancer_sites_filename = ARGV[0]
raise 'Specify random site infos file'  unless random_sites_filename = ARGV[1]

motif_stats = motif_names.map{|motif_name| [motif_name, MotifStat.empty] }.to_h # a: cancer, b: control; positive: core, negative: flank

each_site = ->(filename, &block) {
  # MutatatedSiteInfo.each_site(filename).select{|site| site.disrupted?(fold_change_cutoff: 1) }.each(&block)
  MutatatedSiteInfo.each_site(filename).each(&block)
}

each_site.call(cancer_sites_filename) do |site|
  if site.mutation_in_core?
    motif_stats[site.motif_name].cancer_in_core += 1 # add_a_positive
  else
    motif_stats[site.motif_name].cancer_in_flank += 1 # add_a_negative
  end
end

each_site.call(random_sites_filename) do |site|
  if site.mutation_in_core?
    motif_stats[site.motif_name].random_in_core += 1 # add_b_positive
  else
    motif_stats[site.motif_name].random_in_flank += 1 # add_b_negative
  end
end

pvalue_correction_method = 'fdr'
significance_corrector = PvalueCorrector.new(pvalue_correction_method)


significances = motif_stats.map{|motif_name, infos|
  [motif_name, infos.significance]
}.to_h

significances_corrected = significance_corrector.correct_hash(significances)

puts ['Motif', 'Cancer to random core rate ratio', 'Significances', 'Corrected significances', 'Cancer core rate', 'Random core rate', 'Sites in core (cancer)', 'Sites in flank (cancer)', 'Sites in core (random)', 'Sites in flank (random)'].join("\t")
motif_names.sort.select{|motif_name|
  significances_corrected[motif_name] <= significance_rate_cutoff
}.each do |motif_name|
  infos = motif_stats[motif_name]
  cancer_core_rate = infos.cancer_core_rate
  random_core_rate = infos.random_core_rate
  puts [motif_name,
        infos.cancer_to_random_rate.round(3),
        significances[motif_name],
        significances_corrected[motif_name],
        infos.cancer_core_rate.round(3), 
        infos.random_core_rate.round(3),
        infos.cancer_in_core, infos.cancer_in_flank, infos.random_in_core, infos.random_in_flank
      ].join("\t")
end
