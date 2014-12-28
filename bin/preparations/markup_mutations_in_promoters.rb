require 'interval_notation'

promoter_intervals_by_chromosome = Hash.new{|h,k| h[k] = [] }

File.open('source_data/gene_tss.txt') do |f|
  f.each_line do |line|
    chr, strand, tss = line.chomp.split("\t").last(3)
    tss_pos = tss.to_i
    if strand == '1'
      promoter_intervals_by_chromosome[chr.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - 2000, tss_pos + 500)
    elsif strand == '-1'
      promoter_intervals_by_chromosome[chr.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - 500, tss_pos + 2000)
    end
  end
end

promoters_by_chromosome = promoter_intervals_by_chromosome.map{|chr, intervals|
  [chr, IntervalNotation::Operations.union(intervals)]
}.to_h

File.open('source_data/SUBSTITUTIONS_13Apr2012_snz.txt') do |f|
  puts f.readline
  f.each_line do |line|
    infos = line.chomp.split("\t")
    chr, pos = infos[2], infos[3].to_i

    if promoters_by_chromosome[chr.to_sym].include_position?(pos)
      mut_type = infos[17]
      infos[17] = [mut_type.strip.empty? ? nil : mut_type, 'Promoter'].compact.join(',')
    end
    puts infos.join("\t")
  end
end
