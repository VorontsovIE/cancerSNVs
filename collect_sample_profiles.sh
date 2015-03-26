DIRECTORIES=`ls --directory ./results/motif_statistics/full/any/samples/PD*`
SAMPLES=`basename --multiple $DIRECTORIES`
 

for CONTEXT  in  ${CONTEXTS}; do
  mkdir -p ./results/sample_profiles/${CONTEXT}/disruption/
  mkdir -p ./results/motif_statistics/full/${CONTEXT}/samples/all
  rm -f ./results/sample_profiles/${CONTEXT}/disruption/genome.txt
  rm -f ./results/sample_profiles/${CONTEXT}/disruption/shuffle.txt
  rm -f ./results/sample_profiles/${CONTEXT}/disruption/shuffle_and_genome.txt

  RANDOM_VARIANTS=`find ./results/motif_statistics/full/${CONTEXT}/ -maxdepth 1 -type f | xargs realpath`

  ln -sf ${RANDOM_VARIANTS} ./results/motif_statistics/full/${CONTEXT}/samples/all

  for SAMPLE in $SAMPLES; do

    ls ./results/motif_statistics/full/${CONTEXT}/samples/${SAMPLE}/random_genome_*.csv | ruby sample_profile.rb --sample-name $SAMPLE >> ./results/sample_profiles/${CONTEXT}/disruption/genome.txt
    
    ls ./results/motif_statistics/full/${CONTEXT}/samples/${SAMPLE}/random_shuffle_*.csv | ruby sample_profile.rb --sample-name $SAMPLE >> ./results/sample_profiles/${CONTEXT}/disruption/shuffle.txt

    ls ./results/motif_statistics/full/${CONTEXT}/samples/${SAMPLE}/random_*.csv | ruby sample_profile.rb --sample-name $SAMPLE >> ./results/sample_profiles/${CONTEXT}/disruption/shuffle_and_genome.txt
  done
done
