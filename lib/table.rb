require_relative 'table_header'

class Table
  class Row
    attr_reader :line, :column_headers, :row_name
    def initialize(line, column_headers, row_name)
      @line = line
      @column_headers = column_headers
      @row_name = row_name
    end

    def [](index)
      return @line[index]  if index.is_a?(Fixnum)
      ind = column_headers.header_index(index)
      raise "Element with index #{index} not found"  unless ind
      @line[ind]
    end

    def []=(index, value)
      return @line[index] = value  if index.is_a?(Fixnum)
      ind = column_headers.header_index(index)
      raise "Element with index #{index} not found"  unless ind
      @line[ind] = value
    end

    def size
      column_headers.size
    end
  end


  attr_reader :lines, :header_row, :header_column

  def header_row=(value)
    value = Header.create(value)
    raise "Number of columns(#{number_of_columns}) differs from header row size(#{value.size})"  if value.size != number_of_columns
    @header_row = value
  end

  def header_column=(value)
    value = Header.create(value)
    raise "Number of rows(#{size}) differs from header column size(#{value.size})"  if value.size != size
    @header_column = value
  end

  def initialize(lines, header_row: nil, header_column: nil)
    num_cols = lines.map(&:size).max
    @lines = lines.map{|line| line + Array.new(num_cols - line.size) }

    raise "Number of rows(#{@lines.size}) differs from header column size(#{header_column.size})"  if header_column && header_column.size != @lines.size
    @header_column = header_column ? Header.create(header_column) : EmptyHeader.new(@lines.size)

    raise "Number of columns(#{num_cols}) differs from header row size(#{header_row.size})"  if header_row && num_cols && header_row.size != num_cols
    @header_row = header_row ? Header.create(header_row) :  EmptyHeader.new(num_cols)
  end

  def select_columns(columns)
    rows = each_row.map do |row|
      Row.new(columns.map{|column| row[column] }, header_row.select_columns(columns), row.row_name)
    end
    Table.create_from_rows(rows)
  end

  def self.read(file, with_header_row: false, with_header_column: false)
    if file.respond_to?(:readlines)
      lines = file.readlines
    elsif file.is_a?(String)
      lines = File.readlines(file)
    else
      raise 'file should be either filename(String) or IO-object (or object responding to #readlines)'
    end

    lines = lines.map{|line| line.chomp.split("\t") }
    if with_header_row
      header_row = lines.shift
      header_row.shift  if with_header_column  # both header leave empty cell in the left top corner
    end

    if with_header_column
      header_column = lines.map{|row| row.first }
      lines.each(&:shift)
    end

    Table.new(lines, header_row: header_row, header_column: header_column)
  end

  def self.create_from_rows(rows)
    col_headers = rows.first.column_headers
    raise 'Inconsistent column headers'  unless rows.all?{|row| row.column_headers == col_headers }
    Table.new(rows.map(&:line), header_row: col_headers, header_column: rows.map(&:row_name))
  end

  def size
    lines.size
  end

  def number_of_columns
    lines.first.size
  end

  def output(output_stream, output_header_row: true, output_header_column: true)
    if ! @header_row.empty? && output_header_row
      if ! @header_column.empty? && output_header_column # both row and column
        output_stream.puts("\t" + @header_row.to_s)
      else
        output_stream.puts(@header_row.to_s)
      end
    end
    if ! @header_column.empty? && output_header_column
      @header_column.zip(lines).each do |hdr, line|
        output_stream.puts([hdr, *line].join("\t"))
      end
    else
      lines.each do |line|
        output_stream.puts(line.join("\t"))
      end
    end
    nil
  end

  def add_columns(headers)
    headers.each do |header_sym, header_name|
      add_column([[header_sym, header_name]], :empty)
    end
  end
  def add_column(header, column = :empty)
    if column == :empty
      @header_row += Header.create(header)
      lines.zip( [] ).each do |line, el|
        line << el
      end
    else
      raise 'Column size should be equal to number of lines in the table'  unless column.size == size
      @header_row += Header.create(header)
      lines.zip(column).each do |line, el|
        line << el
      end
    end
  end

  def etc(collection, first_n: 10, separator: ';', etc_sign: '...')
    collection.first(first_n).join(separator) + (collection.size > first_n ? etc_sign : '')
  end
  protected :etc

  def to_s
    header_row_repr = @header_row.etc # header_row ? '[' + etc(header_row) + ']': 'empty'
    header_column_repr = @header_column.etc # header_column ? '[' + etc(header_column) + ']': 'empty'
    "Table<size: #{size} rows x #{number_of_columns} columns; header_row: #{header_row_repr}; header_column: #{header_column_repr}>"
  end

  def inspect
    to_s
  end

  def each_row
    unless block_given?
      iterator = enum_for(:each_row)
      iterator.extend(IteratorWithEnd)
      return iterator
    end
    lines.each_with_index do |line,index|
      yield Row.new(line, header_row, header_column.header_names[index].last)
    end
  end

  def each
    unless block_given?
      iterator = enum_for(:each)
      iterator.extend(IteratorWithEnd)
      return iterator
    end
    lines.each do |line|
      yield line
    end
  end

  def select_rows
    unless block_given?
      iterator = enum_for(:select_rows)
      iterator.extend(IteratorWithEnd)
      return iterator
    end
    # selected_lines = lines.select.each_with_index do |line,index|
    #   yield Row.new(line, header_row, header_column.header_names[index].last)
    # end
    selected_rows = each_row.select do |row|
      yield row
    end
    Table.new(selected_rows.map(&:line), header_row: header_row, header_column: selected_rows.map(&:row_name))
  end

  def sort_rows_by
    unless block_given?
      iterator = enum_for(:sort_rows)
      iterator.extend(IteratorWithEnd)
      return iterator
    end
    # selected_lines = lines.select.each_with_index do |line,index|
    #   yield Row.new(line, header_row, header_column.header_names[index].last)
    # end
    sorted_rows = each_row.sort_by do |row|
      yield row
    end
    Table.new(sorted_rows.map(&:line), header_row: header_row, header_column: sorted_rows.map(&:row_name))
  end

  include Enumerable
end

module IteratorWithEnd
  def eof?
    peek
    false
  rescue StopIteration
    true
  end
end
