motif = {}

snps = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/all_SNPs.txt').drop(1).each do |snp|
  snp_last = snp.split("\t")[-1]
  snp_last_cut = snp_last.to_f
  if snp_last_cut<(1.0/15)
    motif_name = snp.split("\t")[1]

    if motif.has_key?(motif_name)
      motif[motif_name] += 1
    else
      motif[motif_name] = 1
    end
  end
end;


File.open('/home/ilya/iogen/BioSchool-summer2014/SNP/motif_sites_count_less_15.txt', "w") do |f| 
  motif.sort.each do |k, v|
    f.puts "#{k}\t#{v}" 
  end
end

