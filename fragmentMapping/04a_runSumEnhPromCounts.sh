#!/bin/bash
#SBATCH --job-name="sum frags"
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=64G


min_distance=0
max_distance=30000
maxDistkb=`expr $max_distance / 1000`
minDistkb=`expr $min_distance / 1000`
tssUpstream=100
tssDownstream=100
ctrlStrain="366"
expStrains=( "828" "382" "775" "784" "844")


for expStrain in "${expStrains[@]}"
do
  echo "Processing experiment strain: " $expStrain
  enhancerSet="daugherty"
  echo "min distance (converted to kb): " $minDistkb
  echo "max distance (converted to kb): " $maxDistkb
  echo "enhancerSet: " $enhancerSet
  echo "distance upstream of tss: " $tssUpstream
  echo "distance downstream of tss: " $tssDownstream
  Rscript --vanilla ./_04a_sumEnhPromCounts.R $min_distance $max_distance $enhancerSet $tssUpstream $tssDownstream $ctrlStrain $expStrain


  enhancerSet="jaenes"
  echo "min distance (converted to kb): " $minDistkb
  echo "max distance (converted to kb): " $maxDistkb
  echo "enhancerSet: " $enhancerSet
  echo "distance upstream of tss: " $tssUpstream
  echo "distance downstream of tss: " $tssDownstream
  Rscript --vanilla ./_04a_sumEnhPromCounts.R $min_distance $max_distance $enhancerSet $tssUpstream $tssDownstream $ctrlStrain $expStrain

done
