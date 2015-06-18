namespace 'source_data' do
  desc 'Download Wingender TF ontology'
  task :WingenderClassification => 'source_data/TFOntologies/TFClass_human.obo'

  file 'source_data/TFClass_ontologies.zip' do
    sh 'wget', 'http://tfclass.bioinf.med.uni-goettingen.de/suplementary/TFClass_ontologies.zip', '-O', 'source_data/TFClass_ontologies.zip'
  end

  file 'source_data/TFOntologies/TFClass_human.obo' => 'source_data/TFClass_ontologies.zip' do
    sh 'unzip', 'source_data/TFClass_ontologies.zip', '-d', 'source_data/TFOntologies/'
  end
end
