namespace 'source_data' do
  # Nik-Zainal et al. data
  namespace 'YeastApobec' do
    desc 'Unpack somatic mutations data introduced by APOBEC inlined in yeast'
    task :prepare do
      #sh 'wget', '-O', LocalPaths::Secondary::NikZainalSNVsOriginal, 'ftp://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl/SUBSTITUTIONS_13Apr2012_snz.txt'
      sh 'unzip', 'source_data/apobec.zip', '-d', 'source_data/YeastApobec/'
    end
  end
end
