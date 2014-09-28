# Attention! only sequences of length: 25+1+25 are allowed
def cpg_mutation?(sequence)
  (sequence[26] == 'C' && sequence[30] == 'G') || (sequence[26] == 'G' && sequence[24] == 'C')
end

def tpc_mutation?(sequence)
  (sequence[26] == 'C' && sequence[24] == 'T') || (sequence[26] == 'G' && sequence[30] == 'A')
end

def intronic_mutation?(mutation_types)
  mutation_types.split(',').map(&:strip).include?('Intronic')
end

def promoter_mutation?(mutation_types)
  mutation_types.split(',').map(&:strip).include?('Promoter')
end
