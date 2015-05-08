def convert_fasta_to_plain(genome_folder)
  Dir.glob(File.join(genome_folder, '*.fa')).each do |fasta_filename|
    num_entries = 0
    File.open(fasta_filename.ext('.plain'), 'w') do |fw|
      File.open(fasta_filename) do |f|
        f.each_line do |line|
          if line.start_with?('>')
            num_entries += 1
            next
          end
          fw.print(line.strip)
        end
      end
    end
    raise "Error in #{fasta_filename}! More than one entry found!!!"  if num_entries > 1
  end
end

namespace 'source_data' do
  namespace 'genome' do
    desc 'Download and prepare genome assembly (Ensembl, GRCh37.p13)'
    task prepare: LocalPaths::Genome

    file LocalPaths::Genome => [SystemPaths::Genome] do
      rm_rf LocalPaths::Genome # not to put a link into folder in already created folder
      ln_sf SystemPaths::Genome, LocalPaths::Genome
    end

    file SystemPaths::Genome do
      mkdir_p SystemPaths::Genome
      ftp_folder = 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna'
      sh 'wget', '--directory-prefix=#{SystemPaths::Genome}', '#{ftp_folder}/README'
      sh 'wget', '--directory-prefix=#{SystemPaths::Genome}', "#{ftp_folder}/Homo_sapiens.GRCh37.75.dna_sm.chromosome.*.fa.gz"
      sh 'gzip', '-d', *Dir.glob(File.join(SystemPaths::Genome, '*.fa.gz'))
      convert_fasta_to_plain(SystemPaths::Genome)
      rm_f Dir.glob(File.join(SystemPaths::Genome,'*.fa'))
    end
  end
end
