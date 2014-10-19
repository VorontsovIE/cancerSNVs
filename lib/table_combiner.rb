# Column-by-column combiner
class FileColumnCombiner
  attr_reader :mode
  def initialize(mode, filenames)
    @mode = mode
    @files = filenames.map{|filename| File.open(filename) }
    raise "Header mode can be one of:\n :one, :zero, :all"  unless [:one, :zero, :all].include?(mode)
  end

  def each_line
    return enum_for(:each_line)  unless block_given?
    loop do
      yield next_line
    end
  end

  def check_end!
    if @files.all?(&:eof?)
      raise StopIteration
    elsif @files.any?(&:eof?)
      raise "Some of files ended, but not all"
    end
  end

  def next_line
    check_end!
    first_line = @files.first.readline.chomp
    rest_lines = @files.drop(1).map(&:readline).map(&:chomp)
    case mode
    when :all
      [first_line, *rest_lines].join("\t")
    when :one
      header, first_line = first_line.split("\t", 2)
      rest_lines = rest_lines.map{|line| line.split("\t", 2)[1] }
      [header, first_line, *rest_lines].join("\t")
    when :zero
      header, first_line = first_line.split("\t", 2)
      rest_lines = rest_lines.map{|line| line.split("\t", 2)[1] }
      [first_line, *rest_lines].join("\t")
    end
  end
end

# Row-by-row combiner
class FileRowCombiner
  attr_reader :mode
  def initialize(mode, filenames)
    @mode = mode
    @filenames = filenames
    @file_number = 0
    raise "Header mode can be one of:\n :one, :zero, :all"  unless [:one, :zero, :all].include?(mode)
  end

  def yield_first_line_if_necessary(first_row)
    case mode
    when :one
      yield first_row  if @file_number == 1 
    when :zero
    when :all
      yield first_row
    end
  end
  
  def yield_from_stream(stream, &block)
    yield_first_line_if_necessary(stream.readline, &block)
    stream.each_line do |line|
      yield line
    end
  end

  protected :yield_from_stream, :yield_first_line_if_necessary

  def each_line(&block)
    return enum_for(:each_line)  unless block_given?
    @filenames.each do |filename|
      @file_number += 1
      File.open(filename) do |f|
        yield_from_stream(f, &block)
      end
    end
  end
end
