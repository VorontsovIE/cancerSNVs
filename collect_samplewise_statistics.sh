#!/bin/bash

# SITES_FOLDER, FITTING_FOLDER and MOTIF_STATISTICS_FOLDER and so on should be specified in environment variables
cd "$(dirname "$0")"

FITTING_FOLD=1
FOLD_CHANGE=5

# Choose only subset of cancer types
declare -A CANCER_SAMPLES_BY_TYPE
CANCER_SAMPLES_BY_TYPE=( \
                        [ER_plus_ve_HER2_plus_ve]=PD4194a,PD4198a \
                        [ER_plus_ve_HER2_minus_ve]=PD4120a,PD4103a,PD4088a,PD4085a,PD3851a \
                        [ER_minus_ve_HER2_plus_ve]=PD4199a,PD4192a \
                        [triple_negative]=PD4248a,PD4086a,PD4109a \
                        [BRCA1]=PD4107a,PD3890a,PD3905a,PD4005a,PD4006a \
                        [BRCA2]=PD3904a,PD3945a,PD4115a,PD4116a \
                        [PD3851a]=PD3851a \
                        [PD3890a]=PD3890a \
                        [PD3904a]=PD3904a \
                        [PD3905a]=PD3905a \
                        [PD3945a]=PD3945a \
                        [PD4005a]=PD4005a \
                        [PD4006a]=PD4006a \
                        [PD4085a]=PD4085a \
                        [PD4086a]=PD4086a \
                        [PD4088a]=PD4088a \
                        [PD4103a]=PD4103a \
                        [PD4107a]=PD4107a \
                        [PD4109a]=PD4109a \
                        [PD4115a]=PD4115a \
                        [PD4116a]=PD4116a \
                        [PD4120a]=PD4120a \
                        [PD4192a]=PD4192a \
                        [PD4194a]=PD4194a \
                        [PD4198a]=PD4198a \
                        [PD4199a]=PD4199a \
                        [PD4248a]=PD4248a \
                        [ER_plus_ve_HER2_minus_ve_without_PD4120a]=PD4103a,PD4088a,PD4085a,PD3851a \
                        [all_except_PD4120a]=PD3851a,PD3890a,PD3904a,PD3905a,PD3945a,PD4005a,PD4006a,PD4085a,PD4086a,PD4088a,PD4103a,PD4107a,PD4109a,PD4115a,PD4116a,PD4192a,PD4194a,PD4198a,PD4199a,PD4248a \
                      )

for CONTEXT  in  any cpg tpc; do
  mkdir -p  ${SITES_FOLDER}/${CONTEXT}/samples
done

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  ruby bin/preparations/filter_mutations.rb  --cancer-samples ${CANCER_SAMPLES_BY_TYPE[$CANCER_TYPE]}  ${SNV_FOLDER}/SNV_infos_cancer.txt  ${SITES_FOLDER}/any/sites_cancer.txt  >  ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt

  ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                              ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                              --contexts TCN  --mutation-types promoter,intronic  \
                                              >  ${SITES_FOLDER}/tpc/samples/sites_cancer_${CANCER_TYPE}.txt

  ruby bin/preparations/filter_mutations.rb   ${SNV_FOLDER}/SNV_infos_cancer.txt  \
                                              ${SITES_FOLDER}/any/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                              --contexts NCG  --mutation-types promoter,intronic  \
                                              >  ${SITES_FOLDER}/cpg/samples/sites_cancer_${CANCER_TYPE}.txt
done


for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  mkdir -p  $FITTING_FOLDER  ${MOTIF_STATISTICS_FOLDER}/fitting_log/  ${MOTIF_STATISTICS_FOLDER}/slices/samples  ${MOTIF_STATISTICS_FOLDER}/full/samples

  for CONTEXT in any cpg tpc; do
    mkdir -p  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples

    for RANDOM_VARIANT  in  random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
      mkdir -p  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples

      # echo "ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  --fold $FITTING_FOLD  >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log"

      ruby fitting_random_sites.rb  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  \
                                    ${SITES_FOLDER}/${CONTEXT}/sites_${RANDOM_VARIANT}.txt  \
                                    --fold $FITTING_FOLD  \
                                    >  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${RANDOM_VARIANT}.txt  \
                                    2>  ${MOTIF_STATISTICS_FOLDER}/fitting_log/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.log

      ln -f  ${SITES_FOLDER}/${CONTEXT}/samples/sites_cancer_${CANCER_TYPE}.txt  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_cancer.txt
    done
  done
done

for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  for CONTEXT in any cpg tpc; do
    for VARIANT  in  cancer random_genome_13  random_genome_15  random_genome_17  random_shuffle_135  random_shuffle_137  random_shuffle_139; do
      mkdir -p  ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}
      ./generate_all_motif_statistics.sh  ${FITTING_FOLDER}/${CONTEXT}/samples/${CANCER_TYPE}/sites_${VARIANT}.txt  \
                                          ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${VARIANT}  \
                                          0.0005  $FOLD_CHANGE

    done
  done
done



for CANCER_TYPE  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  for CONTEXT in any cpg tpc; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}

    for RANDOM_VARIANT  in  random_shuffle_135  random_shuffle_137  random_shuffle_139  random_genome_13  random_genome_15  random_genome_17; do
      ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/cancer  \
                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}  \
                        ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                        --correction fdr  \
                        >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv

      for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
        for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do
          mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}
          ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv  \
                                  --motif-qualities A,B,C,D  --significance 0.05  \
                                  --${DISRUPTION_OR_EMERGENCE}  --${SUBJECTED_OR_PROTECTED}  \
                                  >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}/${RANDOM_VARIANT}.csv
        done
      done
    done

    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}
    ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_genome_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each_genome.csv
    ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_shuffle_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each_shuffle.csv
    ls ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples/${CANCER_TYPE}/random_*.csv | ruby motif_pvalue_stats.rb  >  ${MOTIF_STATISTICS_FOLDER}/pvalue_statitics/${CONTEXT}/samples/${CANCER_TYPE}/compared_to_each.csv
  done

  for CONTEXT in any cpg tpc; do
    for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
      for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do

        FILTERED_FOLDER=${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}
        COMMON_MOTIFS_FOLDER=${MOTIF_STATISTICS_FOLDER}/common_motifs/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples/${CANCER_TYPE}

        mkdir -p  ${COMMON_MOTIFS_FOLDER}

        ls  ${FILTERED_FOLDER}/random_genome_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_genome.txt
        ls  ${FILTERED_FOLDER}/random_shuffle_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each_shuffle.txt
        ls  ${FILTERED_FOLDER}/random_*.csv | ruby common_motifs.rb  >  ${COMMON_MOTIFS_FOLDER}/compared_to_each.txt
      done
    done
  done
done


for CANCER_TYPE_1  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
  for CONTEXT in any cpg tpc; do
    mkdir -p  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/

    for CANCER_TYPE_2  in  ${!CANCER_SAMPLES_BY_TYPE[@]};  do
      ruby summary.rb   ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE_1}/cancer  \
                        ${MOTIF_STATISTICS_FOLDER}/slices/${CONTEXT}/samples/${CANCER_TYPE_2}/cancer  \
                        ./source_data/motif_names.txt  ./source_data/hocomoco_genes_infos.csv  \
                        --correction fdr  \
                        >  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/${CANCER_TYPE_2}.csv

      for SUBJECTED_OR_PROTECTED  in  subjected  protected; do
        for DISRUPTION_OR_EMERGENCE  in  disruption  emergence; do
          mkdir -p ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}
          ruby filter_summary.rb  ${MOTIF_STATISTICS_FOLDER}/full/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/${CANCER_TYPE_2}.csv  \
                                  --motif-qualities A,B,C,D  --significance 0.05  \
                                  --${DISRUPTION_OR_EMERGENCE}  --${SUBJECTED_OR_PROTECTED}  \
                                  >  ${MOTIF_STATISTICS_FOLDER}/filtered/${SUBJECTED_OR_PROTECTED}/${DISRUPTION_OR_EMERGENCE}/${CONTEXT}/samples_vs_samples/${CANCER_TYPE_1}/${CANCER_TYPE_2}.csv
        done
      done
    done
  done
done
