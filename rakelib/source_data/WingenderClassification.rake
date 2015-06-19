namespace 'source_data' do
  desc 'Download Wingender TF ontology'
  task :WingenderClassification => 'source_data/TFOntologies/TFClass_human.obo'

  file 'source_data/TFClass_ontologies.zip' do
    sh 'wget', 'http://tfclass.bioinf.med.uni-goettingen.de/suplementary/TFClass_ontologies.zip', '-O', 'source_data/TFClass_ontologies.zip'
  end

  file 'source_data/TFOntologies/TFClass_human.obo' => 'source_data/TFClass_ontologies.zip' do
    sh 'unzip', 'source_data/TFClass_ontologies.zip', '-d', 'source_data/TFOntologies/'
  end

  desc 'Download Uniprot ID-AC mapping'
  task :uniprot_id_ac_mapping => 'source_data/human_uniprot.txt'

  file 'source_data/human_uniprot.txt.gz' do
    sh 'wget', 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=organism:%22Homo%20sapiens%20(Human)%20[9606]%22&fil=&format=tab&force=yes&columns=id,entry%20name,genes,protein%20names', '-O', 'source_data/human_uniprot.txt.gz'
  end

  file 'source_data/human_uniprot.txt' => 'source_data/human_uniprot.txt.gz' do
    sh 'gzip', '--decompress', 'source_data/human_uniprot.txt.gz'
  end
end
