raise 'Specify file with sample profiles\' vectors'  unless filename = ARGV[0]

# def sigma(x)
#   Math.exp(x) / (1+Math.exp(x))
# end

def dist(weighted_vect_1, weighted_vect_2)
  vect_1, weight_1 = weighted_vect_1
  vect_2, weight_2 = weighted_vect_2
  raise 'Incompatible sizes'  unless vect_1.size == vect_2.size && vect_1.size == weight_1.size && weight_1.size == weight_2.size
  
  averaged_weights = vect_1.each_index.map{|i| (weight_1[i] + weight_2[i]) / 2.0 }
  vect_1.each_index.map{|i|
    ((vect_1[i] - vect_2[i]) ** 2) * averaged_weights[i]
  }.inject(0.0, &:+) ** 0.5 / averaged_weights.inject(0.0, &:+)**0.5
end

def weighted_distance_matrix(weighted_vector_list)
  weighted_vector_list.map{|weighted_vect_1|
    weighted_vector_list.map{|weighted_vect_2|
      dist(weighted_vect_1, weighted_vect_2)
    }
  }
end

sample_profiles = File.readlines(filename).each_slice(2).map{|weights_line, profile_line|
  sample_name = profile_line.chomp.split("\t", 2)[0]
  profile = profile_line.chomp.split("\t").drop(1).map(&:to_f)
  weights = weights_line.chomp.split("\t").drop(1).map(&:to_f)
  [sample_name, [profile, weights]]
}.to_h

sample_names = sample_profiles.keys
weighted_sample_vectors = sample_names.map{|sample_name| sample_profiles[sample_name] }

matrix = weighted_distance_matrix(weighted_sample_vectors)

puts [nil, *sample_names].join("\t")
sample_names.each_with_index do |sample_name, ind|
  puts [sample_name, *matrix[ind]].join("\t")
end
