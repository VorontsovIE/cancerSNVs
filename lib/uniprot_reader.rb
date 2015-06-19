UniprotEntry = Struct.new(:uniprot_ac, :uniprot_id, :gene_name, :protein_name) do
  def self.from_string(line)
    uniprot_ac, uniprot_id, gene_name, protein_name = line.chomp.split("\t")
    self.new(uniprot_ac, uniprot_id, gene_name, protein_name)
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line.drop(1).map{|line| self.from_string(line) }.each(&block)
    end
  end
end
