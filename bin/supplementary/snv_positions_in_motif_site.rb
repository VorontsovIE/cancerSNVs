$:.unshift File.absolute_path('../../lib', __dir__)
require 'site_info'
require 'measurement_vector'
require 'bioinform'
require 'shellwords'
require 'optparse'
require 'fileutils'

def mean(values)
  values.inject(&:+) / values.size.to_f
end

def stddev(values)
  m = mean(values)
  (values.map{|x| (x - m) ** 2 }.inject(&:+) / (values.size.to_f - 1.0)) ** 0.5
end

motif_collection_folder = './source_data/motif_collection/'

motifs = Dir.glob(File.join(motif_collection_folder, '*.pwm')).map{|fn| Bioinform::MotifModel::PWM.from_file(fn) }
motif_names = motifs.map(&:name).map(&:to_sym)

requested_motifs = motif_names
folder = nil # './results/disruption_position_profile'
normalize = false

OptionParser.new{|opts|
  opts.banner = "Usage: #{opts.program_name} <list of files with sites> [options]\n" +
                "         or\n" +
                "       <list of files with sites> | #{opts.program_name} [options]"
  opts.separator "Options:"
  opts.on('--motifs MOTIFS', 'Specify list of comma-separated motif names') do |value|
    requested_motifs = value.split(',')
  end
  opts.on('--folder FOLDER', 'Specify output folder') do |value|
    folder = value
    FileUtils.mkdir_p(folder)  unless Dir.exist?(folder)
  end
}.parse!(ARGV)

if $stdin.tty?
  sites_filenames = ARGV
else
  sites_filenames = $stdin.read.shellsplit
end



mutation_profiles = sites_filenames.map do |sites_filename|
  mutation_profile_by_motif = {}
  motifs.each do |motif|
    mutation_profile_by_motif[ motif.name.to_sym ] = Array.new(motif.length, 0)
  end

  # MutatatedSiteInfo.each_site(sites_filename).select{|site| site.disrupted?(fold_change_cutoff: 1) }.each do |site|
  MutatatedSiteInfo.each_site(sites_filename).each do |site|
    mutation_profile_by_motif[site.motif_name][site.snv_position_in_site_1_pwm] += 1
  end
  mutation_profile_by_motif
end

def print_info(motif_mutation_profiles, stream)
  normalized_profiles = motif_mutation_profiles.map{|motif_profile|
    motif_profile / motif_profile.inject(0.0, &:+)
  }

  stream.puts mean(normalized_profiles).round(4)
  # stream.puts stddev(normalized_profiles).round(4)
end

requested_motifs.each do |motif|
  motif_mutation_profiles = mutation_profiles.map{|mutation_profile|
    Vector.new(mutation_profile[motif])
  }
  if folder
    File.open(File.join(folder, "#{motif}.txt"), 'w') do |f|
      print_info(motif_mutation_profiles, f)
    end
  else
    $stdout.puts motif
    print_info(motif_mutation_profiles, $stdout)
  end
end
