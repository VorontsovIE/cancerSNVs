snp = File.readlines('/home/ilya/iogen/BioSchool-summer2014/SNP/SUBSTITUTIONS_13Apr2012_snz_promoter_markup.txt').drop(1)
snp_splitted = snp.map{|el| el.split("\t")}
mut = snp_splitted.map{|el| el[17]}.uniq
p mut
