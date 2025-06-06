#!/bin/bash
#SBATCH --job-name="subset frags"
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=32G


min_distance=0
max_distance=30000
maxDistkb=`expr $max_distance / 1000`
minDistkb=`expr $min_distance / 1000`
tssUpstream=100
tssDownstream=100


echo "min kb: " $minDistkb
echo "max kb: " $maxDistkb
Rscript --vanilla ./_03_subsetGI.R $min_distance $max_distance $tssUpstream $tssDownstream
