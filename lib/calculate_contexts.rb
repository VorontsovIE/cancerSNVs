def possible_contexts
  @possible_contexts ||= begin
    results = []
    %w[A C G T].each{|five_prime|
      %w[C T].each{|reference_allele|
        %w[A C G T].each{|three_prime|
          (%w[A C G T] - [reference_allele]).each{|mutated_allele|
            results << "#{five_prime}[#{reference_allele}/#{mutated_allele}]#{three_prime}"
          }
        }
      }
    }
    results
  end
end

def possible_short_contexts
  @possible_short_contexts ||= begin
    results = []
    %w[A C G T].each{|five_prime|
      %w[C T].each{|reference_allele|
        %w[A C G T].each{|three_prime|
          results << "#{five_prime}#{reference_allele}#{three_prime}"
        }
      }
    }
    results
  end
end

def context_counts_in_file(snvs_filename)
  context_counts = Hash.new(0)
  SNVInfo.each_in_file(snvs_filename){|snv|
    expanded_context = SequenceWithSNV.from_string(snv.variant_id.split('@').last).in_pyrimidine_context.to_s
    context_counts[expanded_context] += 1
  }
  context_counts
end

# "short" context is not actually short but unexpanded context
def short_context_counts_in_file(snvs_filename)
  context_counts = Hash.new(0)
  SNVInfo.each_in_file(snvs_filename){|snv|
    short_context = SequenceWithSNV.from_string(snv.variant_id.split('@').last).in_pyrimidine_context.sequence_variant(0)
    context_counts[short_context] += 1
  }
  context_counts
end
