$:.unshift File.absolute_path('../../lib', __dir__)
require 'perfectosape/results'
require 'measurement_vector'
require 'experiment_configuration'
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

motifs = Dir.glob(File.join(MOTIF_COLLECTION_FOLDER, '*.pwm')).map{|fn| Bioinform::MotifModel::PWM.from_file(fn) }
motif_names = motifs.map(&:name).map(&:to_sym)

requested_motifs = motif_names
folder = nil # './results/disruption_position_profile'
normalize = false
expand_region_length = 0

OptionParser.new{|opts|
  opts.banner = "Usage: #{opts.program_name} <list of files with sites> [options]\n" +
                "         or\n" +
                "       <list of files with sites> | #{opts.program_name} [options]"
  opts.separator "Options:"
  opts.on('--motifs MOTIFS', 'Specify list of comma-separated motif names') do |value|
    requested_motifs = value.split(',')
  end
  opts.on('--expand-region LENGTH', 'Specify length of flanks') do |value|
    expand_region_length = Integer(value)
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
    mutation_profile_by_motif[ motif.name.to_sym ] = Array.new(motif.length + 2*expand_region_length, 0)
  end

  # PerfectosAPE::Result.each_in_file(sites_filename).select{|site| site.disrupted?(fold_change_cutoff: 1) }.each do |site|
  PerfectosAPE::Result.each_in_file(sites_filename).each do |site|
    pos = site.snv_position_in_site_1_pwm + expand_region_length # works for expansion procedure defined in such a way that site should ovelap SNV or its expand_region_length vicinity
    raise 'Bad coordinates'  if pos < 0
    mutation_profile_by_motif[site.motif_name][pos] += 1
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
