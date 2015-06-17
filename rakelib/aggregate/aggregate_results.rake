require 'set'

def collect_different_sample_statistics(sample_files, stream:)
  motifs_by_sample = sample_files.map{|header, filename|
    [header, File.readlines(filename).map(&:strip).to_set]
  }
  motif_union = motifs_by_sample.map{|header, motifs| motifs }.inject(Set.new, :|).sort
  stream.puts ['Sample', *motif_union].join("\t")
  motifs_by_sample.each{|header, motifs|
    motif_presence = motif_union.map{|motif| motifs.include?(motif) }
    stream.puts [ header, motif_presence.map{|present| present ? 'X' : ''} ].join("\t")
  }
end

def sample_files(context, protected_or_subjected, disruption_or_emergence)
  ( AlexandrovWholeGenomeCancers.map{|sample|
    [sample, File.join('results/motif_statistics/common/Alexandrov/', sample.to_s, 
                        context.to_s, protected_or_subjected.to_s, disruption_or_emergence.to_s, 'compared_to_each.txt') ]
  } +
  YeastApobecSamples.map{|sample|
    [sample, File.join('results/motif_statistics/common/YeastApobec/', sample.to_s, 
                        context.to_s, protected_or_subjected.to_s, disruption_or_emergence.to_s, 'compared_to_each.txt') ]
  } ).to_h
end

directory 'results/motif_statistics/aggregated/'

desc 'Aggregate common motifs over samples'
task 'aggregate_common_motifs' => ['results/motif_statistics/aggregated/'] do
  [:protected, :subjected].each do |protected_or_subjected|
    [:disruption, :emergence].each do |disruption_or_emergence|
      prep = (protected_or_subjected == :subjected) ? 'to' : 'from'
      File.open("results/motif_statistics/aggregated/#{protected_or_subjected}_#{prep}_#{disruption_or_emergence}_in_any_context.txt", 'w') {|fw|
        collect_different_sample_statistics(sample_files('any', protected_or_subjected, disruption_or_emergence), stream: fw)
      }
    end
  end
end
