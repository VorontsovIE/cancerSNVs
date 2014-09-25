hist = Array.new(16,0)
filename = ARGV[0] # '/home/ilya/iogen/BioSchool-summer2014/SNP/shuffled_mutations/all_res_SNPs_shuffle.txt'
fold = File.open(filename) do |f|
  f.each_line.drop(1).map{|l| l.split("\t")[-1].to_f}
end  

values = fold.map{|el| -Math.log2(el)}
values.each do |value|

  for i in 0..15 do
    if value >= 0.5*i && value < 0.5*(i+1)
      hist[i] += 1
    end
  end  
end

p hist