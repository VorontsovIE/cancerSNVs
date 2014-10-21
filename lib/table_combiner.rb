require_relative 'table'

class TableColumnCombiner
  attr_reader :tables
  def initialize(tables)
    @tables = tables
    @table_iterators = tables.map(&:each)
  end

  def combine
    Table.new(each_line.to_a, header_row: header_row, header_column: header_column)
  end

  def each_line
    return enum_for(:each_line)  unless block_given?
    loop do
      yield next_line
    end
  end

  def header_row
    tables.map{|t| t.header_row }.inject(&:+)
  end

  def header_column
    Table::Header.common_header(tables.map(&:header_column))
  end

  def check_end!
    if @table_iterators.all?(&:eof?)
      raise StopIteration
    elsif @table_iterators.any?(&:eof?)
      raise "Some of files ended, but not all of them"
    end
  end

  def next_line
    check_end!
    @table_iterators.map(&:next).inject(&:+)
  end
end


class TableRowCombiner
  attr_reader :tables
  def initialize(tables)
    @tables = tables
  end

  def combine
    Table.new(tables.flat_map(&:lines), header_row: header_row, header_column: header_column)
  end

  def header_row
    Table::Header.common_header(tables.map(&:header_row))
  end

  def header_column
    tables.map{|t| t.header_column }.inject(&:+)
  end
end
