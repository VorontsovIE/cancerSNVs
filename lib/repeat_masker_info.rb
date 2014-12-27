RepeatMaskerInfo = Struct.new(:score,
  :substitution_percent, :deleted_percent, :inserted_percent,
  :sequence_name, :match_start, :match_end,
  :aftermatch_length, :orientation,
  :repeat_name, :repeat_type,
  :match_margin, :repeat_start_margin, :repeat_end_margin, :id) do

  # http://www.repeatmasker.org/webrepeatmaskerhelp.html#reading
  #
  # 1306 15.6  6.2  0.0 HSU08988  6563  6781  (22462) C  MER7A    DNA/MER2_type    (0)   336   103  12345
  #
  #
  # 1306    = Smith-Waterman score of the match, usually complexity adjusted
  #       The SW scores are not always directly comparable. Sometimes
  #       the complexity adjustment has been turned off, and a variety of
  #       scoring-matrices are used.
  # 15.6    = % substitutions in matching region compared to the consensus
  # 6.2     = % of bases opposite a gap in the query sequence (deleted bp)
  # 0.0     = % of bases opposite a gap in the repeat consensus (inserted bp)
  # HSU08988 = name of query sequence
  # 6563    = starting position of match in query sequence
  # 7714    = ending position of match in query sequence
  # (22462) = no. of bases in query sequence past the ending position of match
  # C       = match is with the Complement of the consensus sequence in the database
  # MER7A   = name of the matching interspersed repeat
  # DNA/MER2_type = the class of the repeat, in this case a DNA transposon 
  #           fossil of the MER2 group (see below for list and references)
  # (0)     = no. of bases in (complement of) the repeat consensus sequence 
  #           prior to beginning of the match (so 0 means that the match extended 
  #           all the way to the end of the repeat consensus sequence)
  # 336    = starting position of match in database sequence (using top-strand numbering)
  # 103    = ending position of match in database sequence
  # 12345  --> id
  #
  def self.from_string(str)
    score,
    substitution_percent, deleted_percent, inserted_percent,
    sequence_name,
    match_start, match_end,
    aftermatch_length,
    orientation,
    repeat_name, repeat_type,
    match_margin, repeat_start_margin, repeat_end_margin, id = str.lstrip.chomp.split(/\s+/)

    RepeatMaskerInfo.new(score.to_i,
      substitution_percent.to_f, deleted_percent.to_f, inserted_percent.to_f,
      sequence_name.to_sym, match_start.to_i, match_end.to_i,
      remove_braces(aftermatch_length).to_i,
      convert_orientation(orientation),
      repeat_name.to_sym, repeat_type.to_sym,
      remove_braces_if_present(match_margin).to_i,
      repeat_start_margin.to_i, repeat_end_margin.to_i, id.to_i)
  end

  def self.remove_braces(str)
    str.match(/\A\((.+)\)\z/)[1]
  end

  def self.remove_braces_if_present(str)
    str.start_with?('(') ? remove_braces(str) : str
  end

  def self.convert_orientation(str)
    if str == '+'
      :+
    elsif str == 'C'
      :-
    else
      raise ArgumentError, 'Complement should be specified by `+` or `C`'
    end
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename).lazy  unless block_given?
    File.open(filename) do |f|
      3.times{ f.readline }
      each_in_stream(f, &block)
    end
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.lazy.map{|line| RepeatMaskerInfo.from_string(line) }.each(&block)    
  end
end
