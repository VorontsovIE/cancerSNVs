$SEQ_FOLDER=$1
$CHUNK_FOLDER=$2
$NUMBER_OF_CORES=$3

mkdir -p  $CHUNK_FOLDER

# Split sequences into equal-size chunks for computation parallelization
for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
  split  --number l/${NUMBER_OF_CORES}  --numeric-suffixes 1  --suffix-length `expr length ${NUMBER_OF_CORES}`  --additional-suffix .txt   ${SEQ_FOLDER}/sequences_${VARIANT}.txt  ${CHUNK_FOLDER}/sequences_${VARIANT}_chunk_
  for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
    mkdir -p ${CHUNK_FOLDER}/core_${SUFFIX}
    mv  ${CHUNK_FOLDER}/sequences_${VARIANT}_chunk_${SUFFIX}.txt  ${CHUNK_FOLDER}/core_${SUFFIX}/sequences_${VARIANT}.txt
  done
done


# Generate scripts for in-parallel PerfectosAPE run. Provide motif collection and PerfectosAPE jar files
echo "" > ${CHUNK_FOLDER}/run_perfectosape_multithread.sh
chmod 755 ${CHUNK_FOLDER}/run_perfectosape_multithread.sh

for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
  ln  ./ape.jar  ${CHUNK_FOLDER}/core_${SUFFIX}/ape.jar
  ln  ./source_data/motif_collection  ${CHUNK_FOLDER}/core_${SUFFIX}/motif_collection
  ln  ./source_data/motif_thresholds  ${CHUNK_FOLDER}/core_${SUFFIX}/motif_thresholds

  # create script from scratch
  echo "" > ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  chmod 755 ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
    echo "java -cp ape.jar ru.autosome.perfectosape.SNPScan  ./motif_collection  ./sequences_${VARIANT}.txt  --fold-change-cutoff 1  --precalc ./motif_thresholds  >  ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  done

  # run on all cores in background
  echo "./core_${SUFFIX}/run_perfectosape.sh &"  >>  ${CHUNK_FOLDER}/run_perfectosape_multithread.sh
done


# Generate scripts for concatenating results of parallel invocations
echo "" > ${CHUNK_FOLDER}/concatenate_results.sh
chmod 755 ${CHUNK_FOLDER}/concatenate_results.sh

for VARIANT  in  cancer  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
  echo "echo \"\" > ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
    echo "cat ./core_${SUFFIX}/sites_${VARIANT}.txt  >>  ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  done
  for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
    echo "rm ./core_${SUFFIX}/sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  done
done
