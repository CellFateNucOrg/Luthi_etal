#!/bin/bash
#SBATCH --job-name="count frag pairs"
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=64G


min_distance=0
max_distance=30000
maxDistkb=`expr $max_distance / 1000`
minDistkb=`expr $min_distance / 1000`

echo "min kb: " $minDistkb
echo "max kb: " $maxDistkb
Rscript --vanilla ./_02_countsToGI.R $min_distance $max_distance
