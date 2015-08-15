$:.unshift './lib'
require 'perfectosape/results_short'
require 'set'

=begin
Dir.glob('recurrent_mutations/SNVs/*').each{|cancer_filename|
  sh 'java', '-Xmx1G', '-cp', 'source_data/ape.jar', 'ru.autosome.perfectosape.SNPScan',
            'recurrent_mutations/motifs_subjectedToDisruption/', cancer_filename,
            '--precalc', 'source_data/motif_thresholds',
            '--compact', {out: cancer_filename.pathmap('recurrent_mutations/sitesSubjectedToDisruption/%n.txt')}, {}
}
=end


# http://manticore.niehs.nih.gov/pavis2/annotate ; assembly hg19 knownGene

class PerfectosAPE::ResultShort
  attr_accessor :sample_type
  def mutation_info
     variant_id.match(/^(?<sample>.+);(?<chromosome>.+):(?<position>\d+)@(?<context>.+)$/)   # 400220;10:119249126@A[T/C]G 
  end
  def chromosome; "chr#{mutation_info[:chromosome]}"; end
  def mutation_position; mutation_info[:position].to_i; end
  def site_position; mutation_position + pos_1; end
  def sample; mutation_info[:sample]; end
  def full_sample_name; "#{sample_type}:#{sample}"; end
  def to_s; "#{full_sample_name};#{chromosome}:#{mutation_position}@#{motif_name}"; end
end

# def site_coordinates_in_file(filename)
#   cancer = File.basename(filename, File.extname(filename))
#   result = PerfectosAPE::ResultShort.each_in_file(filename).to_a
#   result.each{|mut_site| mut_site.sample_type = cancer }
#   result
# end

raise 'Specify check mode (disruption/emergence/all)'  unless check_mode = ARGV[0] # :emergence
raise 'Specify threshold number of mutations on site'  unless threshold_num_mutations_in_site = ARGV[1]
raise 'Specify sites folder'  unless sites_folder = ARGV[2]
check_mode = check_mode.downcase.to_sym
threshold_num_mutations_in_site = threshold_num_mutations_in_site.to_i

case check_mode
when :disruption
  checker = ->(mut_site){ mut_site.disrupted?(fold_change_cutoff: 4) }
when :emergence
  checker = ->(mut_site){ mut_site.emerged?(fold_change_cutoff: 4) }
when :all
  checker = ->(mut_site){ true }
else
  raise 'Unknown mode'
end

mutations_counts_by_site = Hash.new(0)
mutations_by_site = Hash.new{|h,k| h[k] = [] }

Dir.glob(File.join(sites_folder, '*.txt')).each{|filename|
  $stderr.puts(filename)
  PerfectosAPE::ResultShort.each_in_file(filename).lazy.select(&checker).each{|mut_site|
    key = [mut_site.motif_name, mut_site.chromosome, mut_site.site_position]
    mutations_counts_by_site[key] += 1
  }
  mutations_counts_by_site.select!{|k,v| v >= threshold_num_mutations_in_site }
}

Dir.glob(File.join(sites_folder, '*.txt')).each{|filename|
  $stderr.puts(filename)
  cancer = File.basename(filename, File.extname(filename))
  PerfectosAPE::ResultShort.each_in_file(filename).lazy.map{|mut_site|
    mut_site.tap{|x| x.sample_type = cancer }
  }.select(&checker).each{|mut_site|
    key = [mut_site.motif_name, mut_site.chromosome, mut_site.site_position]
    next  unless mutations_counts_by_site[key] >= threshold_num_mutations_in_site
    mutations_by_site[key] << mut_site
  }
  GC.start
}
$stderr.puts 'loaded'

uniq_mutations_by_site = mutations_by_site.map{|site, mutations|
  [site, mutations.uniq(&:full_sample_name)]
}
$stderr.puts 'uniqued'

recurrent_mutations_by_sites = uniq_mutations_by_site.select{|site, mutations|
  mutations.size >= threshold_num_mutations_in_site
}.map{|site, mutations|
  [Set.new([site]),mutations]
}
$stderr.puts 'filtered'

def output_recurrent_mutations(recurrent_mutations_by_sites)
  recurrent_mutations_by_sites.each do |sites, mutations|
    from, to = sites.flat_map{|motif, chr, pos|
      [pos, pos + Bioinform::MotifModel::PWM.from_file("/home/ilya/iogen/hocomoco/#{motif}.pwm").length]
    }.minmax
    chr = sites.first[1]
    puts '=================='
    puts '> ' + [chr, from - 1, to, sites.map{|motif, chr, pos| motif }.join(',')].join("\t")
    puts mutations.map{|mut| "#{mut.full_sample_name};#{mut.chromosome}:#{mut.mutation_position}"}.uniq
  end
end

output_recurrent_mutations(recurrent_mutations_by_sites)
puts '----------------------------------------------'
while true
  $stderr.puts 'glue'
  $stderr.puts recurrent_mutations_by_sites.size
  recurrent_mutations_by_sites.each_with_index{|(sites_1,muts_1), ind_1|
    muts_1.each{|mut|
      recurrent_mutations_by_sites.each_with_index{|(sites_2, muts_2), ind_2|
        next  if (muts_1.map{|mut| [mut.chromosome,mut.mutation_position]} & muts_2.map{|mut| [mut.chromosome,mut.mutation_position]}).empty?
        # p [sites_1, sites_2]
        sites = sites_1 | sites_2
        muts = (muts_1 + muts_2).uniq
        recurrent_mutations_by_sites[ind_1] = [sites, muts]
        recurrent_mutations_by_sites[ind_2] = [sites, muts]
      }
    }
  }
  break  unless recurrent_mutations_by_sites.uniq!
end

output_recurrent_mutations(recurrent_mutations_by_sites)
