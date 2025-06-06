#!/bin/bash
#SBATCH --job-name="count frag pairs"
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-4

source $CONDA_ACTIVATE bioframe

inputFile=( "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20211101_combine_Arima_hicpro_366_3/hicpro/results/hic_results/data/366_3/366_3.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20211129_combine_Arima_hicpro_366_4/hicpro/results/hic_results/data/366_4/366_4.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20211101_combine_Arima_hicpro_828_1/hicpro/results/hic_results/data/828_1/828_1.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20211129_combine_Arima_hicpro_828_2/hicpro/results/hic_results/data/828_2/828_2.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20211101_combine_Arima_hicpro_784_3/hicpro/results/hic_results/data/784_3/784_3.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20211129_combine_Arima_hicpro_784_4/hicpro/results/hic_results/data/784_4/784_4.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20220926_Arima_hicpro_844_1/hicpro/results/hic_results/data/844_1/844_1.allValidPairs" \
    "/mnt/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/20220926_Arima_hicpro_844_2/hicpro/results/hic_results/data/844_2/844_2.allValidPairs" \
)

sampleName=( "366_3" "366_4" "828_1" "828_2" "784_3" "784_4" "844_1" "844_2" )

min_distance=0
max_distance=30000
maxDistkb=`expr $max_distance / 1000`
minDistkb=`expr $min_distance / 1000`

#regionSet=daugherty
#regionSet=jaenes
#regions="./joined${regionSet}EnhancerFragments.bed"

i=`expr $SLURM_ARRAY_TASK_ID - 1`

mkdir -p ./fragCounts
python _01_fragment_pairs2.py -i ${inputFile[$i]} -o ./fragCounts/${sampleName[$i]}_fragment_pair_counts_${minDistkb}-${maxDistkb}kb.txt  \
    -n ${min_distance} -d ${max_distance}

