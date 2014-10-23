# Attention! only sequences of length: 25+1+25 are allowed
def cpg_mutation?(sequence)
  left, before,after, right = sequence.split(/[\[\]\/]/)
  (before == 'C' && right[0] == 'G') || (before == 'G' && left[-1] == 'C')
end

def tpc_mutation?(sequence)
  left, before,after, right = sequence.split(/[\[\]\/]/)
  (before == 'C' && left[-1] == 'T') || (before == 'G' && right[0] == 'A')
end

def intronic_mutation?(mutation_types)
  mutation_types.split(',').map(&:strip).include?('Intronic')
end

def promoter_mutation?(mutation_types)
  mutation_types.split(',').map(&:strip).include?('Promoter')
end


# condition_block: mut_name, sequence --> retain(boolean)
def mutation_names_by_mutation_context(snps_splitted, &condition_block)
  snps_splitted.select(&condition_block).map{|mut_name, sequence| mut_name }.to_set
end

# condition_block: mut_name, mut_type --> retain(boolean)
def mutation_names_by_mutation_type(mutation_types, &condition_block)
  mutation_types.select(&condition_block).map{|mut_name, mut_type| mut_name }.to_set
end
