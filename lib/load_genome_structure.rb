def load_promoters_by_chromosome(filename, length_5_prime: 2000, length_3_prime: 500)
  promoter_intervals_by_chromosome = Hash.new{|h,k| h[k] = [] }

  File.open(filename) do |f|
    f.each_line do |line|
      chromosome, strand, tss = line.chomp.split("\t").last(3)

      chromosome_name = chromosome.to_s.start_with?('CHR') ? chromosome : "chr#{chromosome}"

      tss_pos = tss.to_i
      if strand == '1'
        promoter_intervals_by_chromosome[chromosome_name.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_5_prime, tss_pos + length_3_prime)
      elsif strand == '-1'
        promoter_intervals_by_chromosome[chromosome_name.to_sym] << IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_3_prime, tss_pos + length_5_prime)
      end
    end
  end

  promoter_intervals_by_chromosome.map{|chromosome, intervals|
    [chromosome, IntervalNotation::Operations.union(intervals)]
  }.to_h
end

def load_cage_peaks_by_chromosome(filename, length_5_prime: 2000, length_3_prime: 500)
  CagePeak.each_in_file(filename)
          .group_by(&:chromosome)
          .map{|chromosome, peaks|
            expanded_peak_regions = peaks.map{|peak| peak.region_expanded(length_5_prime: length_5_prime, length_3_prime: length_3_prime) }
            [chromosome, IntervalNotation::Operations.union(expanded_peak_regions)]
          }.to_h
end

def read_repeats_by_chromosome(genome_repeat_masker_folder, ignore_repeat_types: [], expand_length: 0)
  result = {}
  Dir.glob("#{genome_repeat_masker_folder}/*/*.fa.out") do |filename|
    chromosome_name = File.basename(filename, '.fa.out')
    repeats = RepeatMaskerInfo.each_in_file(filename).reject{|info|
      ignore_repeat_types.include?(info.repeat_type)
    }.map{|info|
      IntervalNotation::Syntax::Long.closed_closed(info.match_start - expand_length, info.match_end + expand_length)
    }.to_a # lazy iterator needs to be forced

    result[chromosome_name.to_sym] = IntervalNotation::Operations.union(repeats)
  end
  result
end

def read_coding_exons_by_chromosome(filename)
  EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              coding_regions = exons.map{|exon| exon.coding_part_region }.compact
              chromosome_name = chromosome.to_s.start_with?('chr') ? chromosome : "chr#{chromosome}"
              [chromosome_name.to_sym, IntervalNotation::Operations.union(coding_regions)]
            }.to_h
end

def read_introns_by_chromosome(filename)
  EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              transcript_introns = exons.group_by(&:ensembl_transcript_id).map{|ensembl_transcript_id, exons|
                gene_exons = IntervalNotation::Operations.union(exons.map(&:exon_region))
                gene_exons.covering_interval - gene_exons.complement # introns
              }

              # coding_regions = exons.map{|exon| exon.coding_part_region }.compact
              chromosome_name = chromosome.to_s.start_with?('chr') ? chromosome : "chr#{chromosome}"
              [chromosome_name.to_sym, IntervalNotation::Operations.union(transcript_introns)]
            }.to_h
end
