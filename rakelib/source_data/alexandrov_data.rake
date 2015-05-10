def fix_CRLF(filename)
  file_content = File.read(filename)
  File.write(filename, file_content.gsub("\r","\n"))
end

# Convert Alexandrov et al. supplementary table 4 from xls into csv
def decode_coordinates_of_kataegis(xls_filename, csv_filename)
  require 'roo'
  require 'roo-xls'

  xls = Roo::Excel.new(xls_filename)
  lines = xls.drop(3).map{|line|
    line.map{|el|
      el.is_a?(Numeric) ? el.to_i.to_s : el
    }.join("\t")
  }
  File.write(csv_filename, lines.join("\n"))
end

namespace 'source_data' do
  # Alexandrov et al. data
  namespace 'Alexandrov' do
    desc 'Download somatic mutations data (for several cancer types) from Alexandrov et al.'
    task :prepare do
      rm_rf LocalPaths::Secondary::AlexandrovData
      coordinates_of_kataegis_xls = File.join(LocalPaths::Secondary::AlexandrovData, 'coordinates_of_kataegis.xls')
      sh 'wget', '--recursive', "--directory-prefix=#{LocalPaths::Secondary::AlexandrovData}", 'ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/'
      sh 'wget', '-O', coordinates_of_kataegis_xls, 'http://www.nature.com/nature/journal/v500/n7463/extref/nature12477-s2.xls'

      mv Dir.glob(File.join(LocalPaths::Secondary::AlexandrovData, 'ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/', '*')), LocalPaths::Secondary::AlexandrovData
      fix_CRLF(File.join(LocalPaths::Secondary::AlexandrovData, 'samples_summary.txt'))
      decode_coordinates_of_kataegis(coordinates_of_kataegis_xls, LocalPaths::Secondary::CoordinatesOfKataegis)
      
      rm_rf File.join(LocalPaths::Secondary::AlexandrovData, '/ftp.sanger.ac.uk/')
      rm_f coordinates_of_kataegis_xls
    end
  end
end
