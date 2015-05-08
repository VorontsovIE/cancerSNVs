namespace 'source_data' do
  # Nik-Zainal et al. data
  namespace 'NikZainal' do
    desc 'Download somatic mutations data (for breast cancer only) from Nik-Zainal et al.'
    task :prepare do
      sh 'wget', '-O', LocalPaths::Secondary::NikZainalSNVs, 'ftp://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl/SUBSTITUTIONS_13Apr2012_snz.txt'
    end
  end
end
