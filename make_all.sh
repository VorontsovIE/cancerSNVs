#!/bin/bash

rm -r results/intermediate/site_subsets

mkdir -p results/intermediate/site_subsets


SNV_INFOS = ./source_data/SNV_infos_regulatory.txt
SITES_CANCER = ./results/intermediate/site_subsets/cancer_any.txt
SITES_RANDOM = ./results/intermediate/site_subsets/random_any.txt

# ##  Filtering was already done at SNV-infos stage
# ruby bin/preparations/filter_mutations.rb $SNV_INFOS ./source_data/sites_cancer.txt  --mutation-types promoter,intronic  >  $SITES_CANCER
# ruby bin/preparations/filter_mutations.rb $SNV_INFOS ./source_data/sites_random.txt  --mutation-types promoter,intronic  >  $SITES_RANDOM

ln ./source_data/sites_cancer.txt  $SITES_CANCER
ln ./source_data/sites_random.txt  $SITES_RANDOM

./make_all_for_context.sh  any  $SITES_CANCER  $SITES_RANDOM

ruby bin/preparations/filter_mutations.rb $SNV_INFOS  $SITES_CANCER  --contexts TCN  --mutation-types promoter,intronic  >  ./results/intermediate/site_subsets/cancer_tpc.txt
ruby bin/preparations/filter_mutations.rb $SNV_INFOS  $SITES_RANDOM  --contexts TCN  --mutation-types promoter,intronic  >  ./results/intermediate/site_subsets/random_tpc.txt
./make_all_for_context.sh  tpc  ./results/intermediate/site_subsets/cancer_tpc.txt  ./results/intermediate/site_subsets/random_tpc.txt

ruby bin/preparations/filter_mutations.rb $SNV_INFOS  $SITES_CANCER  --contexts NCG  --mutation-types promoter,intronic  >  ./results/intermediate/site_subsets/cancer_cpg.txt
ruby bin/preparations/filter_mutations.rb $SNV_INFOS  $SITES_RANDOM  --contexts NCG  --mutation-types promoter,intronic  >  ./results/intermediate/site_subsets/random_cpg.txt
./make_all_for_context.sh  cpg  ./results/intermediate/site_subsets/cancer_cpg.txt  ./results/intermediate/site_subsets/random_cpg.txt
