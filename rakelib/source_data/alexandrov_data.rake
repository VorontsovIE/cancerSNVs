def fix_CRLF(filename)
  file_content = File.read(filename)
  File.write(filename, file_content.gsub("\r","\n"))
end

namespace 'source_data' do
  # Alexandrov et al. data
  namespace 'Alexandrov' do
    desc 'Download somatic mutations data (for several cancer types) from Alexandrov et al.'
    task :prepare do
      rm_rf LocalPaths::Secondary::AlexandrovData
      sh 'wget', '--recursive', "--directory-prefix=#{LocalPaths::Secondary::AlexandrovData}", 'ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/'

      mv Dir.glob(File.join(LocalPaths::Secondary::AlexandrovData, 'ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/', '*')), LocalPaths::Secondary::AlexandrovData
      fix_CRLF(File.join(LocalPaths::Secondary::AlexandrovData, 'samples_summary.txt'))
      
      rm_rf File.join(LocalPaths::Secondary::AlexandrovData, '/ftp.sanger.ac.uk/')
    end
  end
end
