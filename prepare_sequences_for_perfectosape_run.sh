#!/bin/bash

# NUMBER_OF_CORES, SEQ_FOLDER and CHUNK_FOLDER should be specified in environment variables

cd "$(dirname "$0")"

mkdir -p  $CHUNK_FOLDER

# Split sequences into equal-size chunks for computation parallelization
for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
  split  --number=l/${NUMBER_OF_CORES}  --numeric-suffixes=1  --suffix-length=`expr length ${NUMBER_OF_CORES}`  --additional-suffix=.txt   ${SEQ_FOLDER}/sequences_${VARIANT}.txt  ${CHUNK_FOLDER}/sequences_${VARIANT}_chunk_
  for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
    mkdir -p ${CHUNK_FOLDER}/core_${SUFFIX}
    mv  ${CHUNK_FOLDER}/sequences_${VARIANT}_chunk_${SUFFIX}.txt  ${CHUNK_FOLDER}/core_${SUFFIX}/sequences_${VARIANT}.txt
  done
done


# Generate scripts for in-parallel PerfectosAPE run. Provide motif collection and PerfectosAPE jar files
echo '#!/bin/bash' > ${CHUNK_FOLDER}/run_perfectosape_multithread.sh
chmod 755 ${CHUNK_FOLDER}/run_perfectosape_multithread.sh
echo 'cd "$(dirname "$0")"'  >>  ${CHUNK_FOLDER}/run_perfectosape_multithread.sh

for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
  ln -f  ./ape.jar  ${CHUNK_FOLDER}/core_${SUFFIX}/ape.jar
  ln -f  ./source_data/motif_collection  ${CHUNK_FOLDER}/core_${SUFFIX}/motif_collection
  ln -f  ./source_data/motif_thresholds  ${CHUNK_FOLDER}/core_${SUFFIX}/motif_thresholds

  # create script from scratch
  echo '#!/bin/bash' > ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  chmod 755 ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  echo 'cd "$(dirname "$0")"'  >>  ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
    echo "java -cp ape.jar ru.autosome.perfectosape.SNPScan  ./motif_collection  ./sequences_${VARIANT}.txt  --fold-change-cutoff 1  --precalc ./motif_thresholds  >  ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/core_${SUFFIX}/run_perfectosape.sh
  done

  # run on all cores in background
  echo "./core_${SUFFIX}/run_perfectosape.sh &"  >>  ${CHUNK_FOLDER}/run_perfectosape_multithread.sh
done
echo "wait"  >>  ${CHUNK_FOLDER}/run_perfectosape_multithread.sh
echo "./concatenate_results.sh"  >>  ${CHUNK_FOLDER}/run_perfectosape_multithread.sh


# Generate scripts for concatenating results of parallel invocations
echo '#!/bin/bash' > ${CHUNK_FOLDER}/concatenate_results.sh
chmod 755 ${CHUNK_FOLDER}/concatenate_results.sh
echo 'cd "$(dirname "$0")"'  >>  ${CHUNK_FOLDER}/concatenate_results.sh

for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
  echo "cat /dev/null > ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  
  # header
  CORE_INDEX=`seq --equal-width  1 ${NUMBER_OF_CORES} | head -1`
  echo "grep -P ^# ./core_${CORE_INDEX}/sites_${VARIANT}.txt  >>  ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  
  for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
    echo "grep --invert-match -P ^# ./core_${SUFFIX}/sites_${VARIANT}.txt  >>  ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh

    # echo "cat ./core_${SUFFIX}/sites_${VARIANT}.txt  >>  ./sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  done
  for SUFFIX in `seq --equal-width  1 ${NUMBER_OF_CORES}`; do
    echo "rm ./core_${SUFFIX}/sites_${VARIANT}.txt"  >>  ${CHUNK_FOLDER}/concatenate_results.sh
  done
done
