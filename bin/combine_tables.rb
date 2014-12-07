$:.unshift File.absolute_path('lib', __dir__)
require 'file_combiner'
require 'optparse'

mode = :all
rows_or_columns = :columns
OptionParser.new do |opts|
  opts.on('--header [MODE]',  "Has header column (sic!); Mode is one/zero/all. If MODE not specifed, `one` is supposed.\n" +
                              "If option not specified, tool works like `all` is supplied") {|value|
    mode = value ? value.to_sym : :one
    raise "Header mode can be one of: 'one', 'zero', 'all'."  unless [:one, :zero, :all].include?(mode)
  }
  opts.on('--rows', "Glue rows, not columns") {
    rows_or_columns = :rows
  }
end.parse!(ARGV)

files = ARGV
raise "Specify at least one file"  unless files.size >= 1

case rows_or_columns
when :columns
  combiner = FileColumnCombiner.new(mode, files)
when :rows
  combiner = FileRowCombiner.new(mode, files)
end

combiner.each_line do |line|
  puts line
end
