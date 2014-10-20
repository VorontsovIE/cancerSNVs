class Table
  attr_reader :lines, :header_row, :header_column

  def initialize(lines, header_row: nil, header_column: nil)
    @lines = lines
    @header_row = header_row
    @header_column = header_column
  end

  def self.read(file, with_header_row: false, with_header_column: false)
    if file.respond_to?(:readlines)
      lines = file.readlines
    elsif lines.is_a? String
      lines = File.readlines(file)
    else
      raise 'file should be either filename(String) or IO-object (or object responding to #readlines)'
    end

    lines = lines.map{|line| line.chomp.split("\t") }
    if with_header_row
      header_row = lines.shift
    end

    if with_header_column
      header_column = lines.map{|row| row.first }
      lines.each(&:shift)
    end

    Table.new(lines, header_row: header_row, header_column: header_column)
  end

  def size
    lines.size
  end

  def output(output_stream)
    output_stream.puts(header_row.join("\t"))  if header_row
    if header_column
      header_column.zip(lines).each do |hdr, line|
        output_stream.puts([hdr, *line].join("\t"))
      end
    else
      lines.each do |line|
        output_stream.puts(line.join("\t"))
      end
    end
  end

  def add_column(header, column)
    raise 'Column size should be equal to number of lines in the table'  unless column.size == size
    @header_row << header  if header_row
    column.each_with_index do |el, row_index|
      @lines[row_index] << el
    end
  end

  # accepts block
  def calculate_column(header)
    @header_row << header  if header_row
    @lines.size.times do |index|
      @lines[index] << (yield @lines[index])
    end
  end

  def each
    return enum_for(:each)  unless block_given?
    lines.each do |line|
      yield line
    end
  end

  include Enumerable
end
