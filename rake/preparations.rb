require 'tmpdir'

namespace 'preparations' do

  namespace 'extractSNVs' do
    desc 'Convert Alexandrov\'s mutations to a common format. Convert format, filter regulatory only SNVs, remove duplicates.'
    task 'AlexandrovEtAl' do
      mkdir_p LocalPaths::Secondary::SNVs
      Dir.mktmpdir('Alexandrov_SNV_infos_temp') do |tmp_dir|
        ruby File.join(LocalPaths::Root, 'bin/preparations/load_cancer_mutations_sequences.rb'), \
             tmp_dir, \
             '--promoter-upstream', '5000', \
             '--promoter-downstream', '500', \
             '--kataegis-expansion', '1000'

        # remove duplicates and put into separate folders
        Dir.glob(File.join(tmp_dir, '*.txt')) do |fn|
          cancer_type = File.basename(fn, '.txt')
          folder = File.join(LocalPaths::Secondary::SNVs, cancer_type)
          mkdir_p folder
          ruby File.join(LocalPaths::Root, 'bin/preparations/filter_snv_infos.rb'), \
               fn, \
               { out: File.join(folder, "#{cancer_type}.txt") }, \
               {} # FileUtils#ruby options :noop/:verbose
        end
      end
    end

    desc 'Convert Nik-Zainal\'s mutations to a common format. Convert format, filter regulatory only SNVs, remove duplicates.'
    task 'NikZainalEtAl' do
      Dir.mktmpdir('NikZainal_SNV_infos_temp') do |tmp_dir|
        ruby 'bin/preparations/convert_breast_cancer_snvs_to_snv_infos.rb', \
             './source_data/SNV_infos_original.txt', \
             { out: File.join(tmp_dir, 'SNV_infos_original.txt') }, \
             {} # FileUtils#ruby options :noop/:verbose

        ruby 'bin/preparations/snv_markup.rb', \
             File.join(tmp_dir, 'SNV_infos_original.txt'), \
             { out: File.join(tmp_dir, 'SNV_infos_marked_up.txt') }, \
             {} # FileUtils#ruby options :noop/:verbose

        ruby 'bin/preparations/filter_snv_infos.rb', \
             File.join(tmp_dir, 'SNV_infos_marked_up.txt'), \
             { out: File.join(LocalPaths::Secondary::SNVs, 'SNV_infos_cancer.txt') }, \
             {} # FileUtils#ruby options :noop/:verbose
      end
    end
  end

end


# mkdir_p LocalPaths::Secondary::Sequences
# mkdir_p LocalPaths::Secondary::Chunks

