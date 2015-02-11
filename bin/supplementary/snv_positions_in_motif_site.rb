$:.unshift File.absolute_path('../../lib', __dir__)
require 'site_info'
require 'bioinform'
require 'shellwords'
require 'optparse'
require 'fileutils'

class Vector
  attr_reader :values
  def initialize(values)
    @values = values
  end

  def [](ind)
    @values[ind]
  end

  def dim
    values.size
  end

  def +(other)
    raise 'Incompatible dimensionality'  unless dim == other.dim
    Vector.new( dim.times.map{|i| self[i] + other[i]} )
  end

  def -(other)
    raise 'Incompatible dimensionality'  unless dim == other.dim
    Vector.new( dim.times.map{|i| self[i] - other[i]} )
  end

  def *(other)
    if other.is_a?(Vector)
      raise 'Incompatible dimensionality'  unless dim == other.dim
      Vector.new( dim.times.map{|i| self[i] * other[i] } )
    elsif other.is_a?(Numeric)
      Vector.new( @values.map{|el| el * other } )
    else
      raise
    end
  end

  def /(other)
    if other.is_a?(Vector)
      raise 'Incompatible dimensionality'  unless dim == other.dim
      Vector.new( dim.times.map{|i| self[i] / other[i] } )
    elsif other.is_a?(Numeric)
      Vector.new( @values.map{|el| el / other } )
    else
      raise
    end
  end

  def **(k)
    Vector.new( @values.map{|el| el ** k } )
  end

  def coerce(other)
    if other.is_a?(Numeric)
      [self, other]
    else
      raise
    end
  end

  def each(&block)
    @values.each(&block)
  end

  def to_s
    @values.join("\t")
  end

  def inspect
    '<' + @values.join(',') + '>'
  end

  def round(rate = 0)
    Vector.new(@values.map{|el| el.round(rate) })
  end

  include Enumerable
end

def mean(values)
  values.inject(&:+) / values.size.to_f
end

def stddev(values)
  m = mean(values)
  (values.map{|x| (x - m) ** 2 }.inject(&:+) / (values.size.to_f - 1.0)) ** 0.5
end

motif_collection_folder = '/home/ilya/iogen/hocomoco/'

motifs = Dir.glob(File.join(motif_collection_folder, '*.pwm')).map{|fn| Bioinform::MotifModel::PWM.from_file(fn) }
motif_names = motifs.map(&:name)

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
    mutation_profile_by_motif[ motif.name ] = Array.new(motif.length, 0)
  end

  MutatatedSiteInfo.each_site(sites_filename).each do |site|
    mutation_profile_by_motif[site.motif_name][site.snv_position_in_site_1_pwm] += 1
  end
  mutation_profile_by_motif
end

def print_info(motif_mutation_profiles, stream_main, stream_additional)
  motif_mutation_profiles.each{|motif_profile|  stream_additional.puts motif_profile }
  stream_additional.puts '========================'


  normalized_profiles = motif_mutation_profiles.map{|motif_profile|
    motif_profile / motif_profile.inject(0.0, &:+)
  }

  normalized_profiles.each{|motif_profile|  stream_additional.puts motif_profile.round(4) }
  stream_additional.puts '========================'

  stream_main.puts mean(normalized_profiles).round(4)
  stream_main.puts stddev(normalized_profiles).round(4)
  stream_additional.puts (100 * stddev(normalized_profiles) / mean(normalized_profiles)).round(2)
end

requested_motifs.each do |motif|
  motif_mutation_profiles = mutation_profiles.map{|mutation_profile|
    Vector.new(mutation_profile[motif])
  }
  if folder
    $stderr.puts motif
    File.open(File.join(folder, "#{motif}.txt"), 'w') do |f|
      print_info(motif_mutation_profiles, f, $stderr)
    end
  else
    puts motif
    print_info(motif_mutation_profiles, $stdout, $stderr)
  end
end
