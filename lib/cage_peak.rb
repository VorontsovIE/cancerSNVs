require 'interval_notation'

CagePeak = Struct.new(:annotation, :short_description, :description, :association_with_transcript,
                  :entrezgene_id, :hgnc_ids, :uniprot_id) do

  attr_reader :chromosome, :strand, :pos_start, :pos_end

  def chromosome
    @chromosome ||= annotation.split(':').first.to_sym
  end

  def strand
    @strand ||= begin
      chromosome, name = annotation.split(':')
      result = name.split(',').last.to_sym
      raise "Unknown strand #{result}"  unless [:+, :-].include?(result)
      result
    end
  end

  def pos_start
    @pos_start ||= begin
      chromosome, name = annotation.split(':')
      name, strand = name.split(',')
      name.split(/\.\.|-/).first.to_i
    end
  end

  def pos_end
    @pos_end ||= begin
      chromosome, name = annotation.split(':')
      name, strand = name.split(',')
      name.split(/\.\.|-/).last.to_i
    end
  end

  def region
    @region ||= IntervalNotation::Syntax::Long.closed_open(pos_start, pos_end)
  end

  def region_expanded(length_5_prime: 2000, length_3_prime: 500)
    case strand
    when :+
      IntervalNotation::Syntax::Long.closed_open(pos_start - length_5_prime, pos_end + length_3_prime)
    when :-
      IntervalNotation::Syntax::Long.closed_open(pos_start - length_3_prime, pos_end + length_5_prime)
    end
  end

# (!) human FANTOM
  def self.new_by_infos(infos)
    annotation, short_description, description, association_with_transcript, entrezgene, hgnc, uniprot_id, *_tpms = infos.chomp.split("\t")
    hgnc_ids = hgnc.split(',')
    entrezgene_ids = entrezgene.split(',')
    CagePeak.new(annotation, short_description, description, association_with_transcript, entrezgene_ids, hgnc_ids, uniprot_id)
  end

# # (!) mouse FANTOM
#   def self.new_by_infos(infos)
#     annotation, short_description, description, association_with_transcript, entrezgene, uniprot_id, tpm = infos.chomp.split("\t")
#     tpm = tpm.to_f
#     # hgnc_ids = hgnc.split(',').map{|hgnc_id| hgnc_from_string(hgnc_id)}
#     entrezgene_ids = entrezgene.split(',').map{|entrezgene_id| entrezgene_from_string(entrezgene_id)}
#     Peak.new(annotation, short_description, description, association_with_transcript, entrezgene_ids, [], uniprot_id, tpm)
#   end

  def self.each_in_stream(stream, &block)
    return enum_for(:each_in_stream, stream).lazy  unless block_given?
    stream.each_line.lazy.select{|line|
      line.start_with?('chr')
    }.map{|line|
      CagePeak.new_by_infos(line)
    }.each(&block)
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename).lazy  unless block_given?
    File.open(filename) do |f|
      CagePeak.each_in_stream(f, &block)
    end
  end
end
