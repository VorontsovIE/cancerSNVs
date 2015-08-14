$:.unshift './lib'
require 'perfectosape/results_short'

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
end

def site_coordinates_in_file(filename)
  cancer = File.basename(filename, File.extname(filename))
  result = PerfectosAPE::ResultShort.each_in_file(filename).to_a
  result.each{|mut_site| mut_site.sample_type = cancer }
  result
end

threshold_num_mutations_in_site = (ARGV[0] || 4).to_i

mutations_by_site = Hash.new{|h,k| h[k] = [] }
Dir.glob('recurrent_mutations/sites/*.txt')
  .flat_map{|fn| site_coordinates_in_file(fn) }
  .select{|mut_site| mut_site.disrupted? }
  .each{|mut_site|
    key = [mut_site.motif_name, mut_site.chromosome, mut_site.site_position]
    mutations_by_site[key] << mut_site
  }

mutations_by_site_uniq = mutations_by_site.map{|site, mutations|
  [site, mutations.uniq(&:full_sample_name)]
}.to_h

mutations_by_site_uniq.select{|site, mutations|
  mutations.size >= threshold_num_mutations_in_site
}.flat_map{|site, mutations|
  mutations.map{|mut|
    [mut.chromosome, mut.mutation_position]
  }.uniq
}.map{|chr, pos|
  [chr, pos - 1, pos]
}.uniq.each{|bed_line|
  puts bed_line.join("\t")
}
