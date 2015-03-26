#!/bin/bash

# SNV_FOLDER, CHUNK_FOLDER and SITES_FOLDER should be specified in environment variables
# CONTEXT_DECLARATIONS should be declared as an associative array
cd "$(dirname "$0")"

declare -A CONTEXT_DECLARATIONS # All except any
CONTEXT_DECLARATIONS=( \
                        [tpc]="--contexts TCN  --mutation-types promoter,intronic" \
                        [cpg]="--contexts NCG  --mutation-types promoter,intronic" \
                        [non_tpc]="--invert-context-request  --contexts TCN  --mutation-types promoter,intronic" \
                      )

for CONTEXT in ${CONTEXTS}; do
  mkdir -p ${SITES_FOLDER}/${CONTEXT}
done

for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do
  ln -f  ${CHUNK_FOLDER}/sites_${VARIANT}.txt  ${SITES_FOLDER}/any/sites_${VARIANT}.txt
done

# generate subsets of sites in TpC/CpG mutation contexts
for CONTEXT in ${CONTEXTS}; do
  if [[ "$CONTEXT" != "any" ]]; then
    for VARIANT  in  cancer  ${RANDOM_VARIANTS}; do

      ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_${VARIANT}.txt  \
                                                  ${SITES_FOLDER}/any/sites_${VARIANT}.txt  \
                                                  ${CONTEXT_DECLARATIONS[$CONTEXT]} \
                                                  >  ${SITES_FOLDER}/${CONTEXT}/sites_${VARIANT}.txt
    done
  fi
done
