hocomoco_archive_file = File.join(LocalPaths::Root, 'source_data/hocomoco_v9_pwm.tar.gz')
hocomoco_pcm_archive_file = File.join(LocalPaths::Root, 'source_data/hocomoco_v9_pcm.tar.gz')

namespace 'source_data' do
  namespace 'motif_collection' do
    desc 'Download and prepare TFBS model collection'
    task prepare: [:motifs, :calculate_thresholds]

    task motifs: [LocalPaths::MotifCollection, LocalPaths::MotifCollectionPCM, LocalPaths::Secondary::MotifNames, LocalPaths::Secondary::HocomocoUniprots, LocalPaths::Secondary::MotifQualities]
    file LocalPaths::MotifCollection => [SystemPaths::MotifCollection] do
      rm_rf LocalPaths::MotifCollection # not to put a link into folder in already created folder
      ln_sf SystemPaths::MotifCollection, LocalPaths::MotifCollection
    end
    file LocalPaths::MotifCollectionPCM => [SystemPaths::MotifCollectionPCM] do
      rm_rf LocalPaths::MotifCollectionPCM # not to put a link into folder in already created folder
      ln_sf SystemPaths::MotifCollectionPCM, LocalPaths::MotifCollectionPCM
    end
    file LocalPaths::Secondary::MotifNames => [LocalPaths::MotifCollection] do
      motif_names = Dir.glob(File.join(LocalPaths::MotifCollection, '*.pwm')).map{|fn| File.basename(fn, '.pwm') }.sort
      File.write(LocalPaths::Secondary::MotifNames, motif_names.join("\n"))
    end

    file LocalPaths::Secondary::MotifQualities do
      sh 'wget', '-O', LocalPaths::Secondary::MotifQualities, 'http://opera.autosome.ru/downloads/motif_collections/hocomoco_qualities.tsv'
    end

    file LocalPaths::Secondary::HocomocoUniprots do
      sh 'wget', '-O', LocalPaths::Secondary::HocomocoUniprots, 'http://opera.autosome.ru/downloads/motif_collections/HOCOMOCOv9_motifs2uniprot.txt'
    end

    file SystemPaths::MotifCollection do
      sh 'wget', '-O', hocomoco_archive_file, 'http://opera.autosome.ru/downloads/motif_collections/hocomoco_v9_pwm.tar.gz'
      mkdir_p  SystemPaths::MotifCollection
      sh 'tar', "--directory=#{SystemPaths::MotifCollection}", '-zxf', hocomoco_archive_file
      rm_f  hocomoco_archive_file
    end

    file SystemPaths::MotifCollectionPCM do
      sh 'wget', '-O', hocomoco_pcm_archive_file, 'http://opera.autosome.ru/downloads/motif_collections/hocomoco_v9_pcm.tar.gz'
      mkdir_p  SystemPaths::MotifCollectionPCM
      sh 'tar', "--directory=#{SystemPaths::MotifCollectionPCM}", '-zxf', hocomoco_pcm_archive_file
      rm_f  hocomoco_pcm_archive_file
    end

    task calculate_thresholds: [LocalPaths::MotifThresholds]
    file LocalPaths::MotifThresholds => SystemPaths::MotifThresholds do
      rm_rf LocalPaths::MotifThresholds # not to put a link into folder in already created folder
      ln_sf SystemPaths::MotifThresholds, LocalPaths::MotifThresholds
    end

    file SystemPaths::MotifThresholds => [:motifs, SystemPaths::PerfectosAPE] do
      sh 'java', '-cp', SystemPaths::PerfectosAPE, 'ru.autosome.ape.PrecalculateThresholds', SystemPaths::MotifCollection, SystemPaths::MotifThresholds, '--silent'
    end
  end
end
