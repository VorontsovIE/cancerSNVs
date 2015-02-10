$:.unshift File.absolute_path('../../lib', __dir__)
require 'site_info'
require 'bioinform'

raise 'Specify file with sites'  unless sites_filename = ARGV[0]

motif_collection_folder = '/home/ilya/iogen/hocomoco/'

mutation_profile_by_motif = {}
Dir.glob(File.join(motif_collection_folder, '*.pwm')).each do |fn|
  mutation_profile_by_motif[ File.basename(fn, File.extname(fn)) ] = Array.new(Bioinform::MotifModel::PWM.from_file(fn).length, 0)
end


MutatatedSiteInfo.each_site(sites_filename).each{|site|
  # mutation_profile_by_motif[site.motif_name][site.snv_position_in_site_1_pwm] ||= 0
  mutation_profile_by_motif[site.motif_name][site.snv_position_in_site_1_pwm] += 1
}

mutation_profile_by_motif.sort.each do |motif, snv_positions|
  puts [motif, snv_positions.inspect].join("\t")
end
