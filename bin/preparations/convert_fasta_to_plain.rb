dirname = ARGV.shift

Dir.glob(File.join(dirname, '*.fa')).sort.each do |filename|
  basename = File.basename(filename, File.extname(filename))
  dirname = File.dirname(filename)
  filename_plain = File.join(dirname, "#{basename}.plain")
  $stderr.print "#{filename} --> #{filename_plain} started"
  num_entries = 0
  File.open(filename_plain, 'w') do |fw|
    File.open(filename) do |f|
      f.each_line do |line|
        if line.start_with?('>')
          num_entries += 1
          next
        end
        fw.print line.strip
      end
    end
  end
  $stderr.print " --> complete\n"
  $stderr.puts "Error in #{filename}! More than one entry found!!!"  if num_entries > 1
end
