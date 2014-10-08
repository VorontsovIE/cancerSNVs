require 'optparse'

class FileColumnCombiner
  attr_reader :mode, :files
  def initialize(mode, filenames)
    @mode = mode
    @files = filenames.map{|filename| File.open(filename) }
    raise "Header mode can be one of:\n :one, :zero, :all"  unless [:one, :zero, :all].include?(mode)
  end

  def each_line
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

mode = :all
OptionParser.new do |opts|
  opts.on('--header [MODE]',  "Has header column (sic!); Mode is one/zero/all. If MODE not specifed, `one` is supposed.\n" +
                              "If option not specified, tool works like `all` is supplied") {|value|
    mode = value ? value.to_sym : :one
    raise "Header mode can be one of: 'one', 'zero', 'all'."  unless [:one, :zero, :all].include?(mode)
  }
end.parse!(ARGV)

files = ARGV
raise "Specify at least one file"  unless files.size >= 1

combiner = FileColumnCombiner.new(mode, files)
combiner.each_line do |line|
  puts line
end
