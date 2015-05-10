namespace 'source_data' do
  desc 'Download PerfectosAPE jar-file'
  task :perfectosape => LocalPaths::PerfectosAPE

  file LocalPaths::PerfectosAPE => [SystemPaths::PerfectosAPE] do
    ln_sf SystemPaths::PerfectosAPE, LocalPaths::PerfectosAPE
  end

  file SystemPaths::PerfectosAPE do
    sh 'wget', 'http://opera.autosome.ru/downloads/ape.jar', '-O', SystemPaths::PerfectosAPE
  end
end
