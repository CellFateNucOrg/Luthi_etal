#! /bin/bash
micromamba activate deeptools

COH1=/Volumes/external.data/MeisterLab/bisiaka/supercoilingdata/supercoiling2/fountainspaperfiles/bigwigs/COH-1modencode_mean.bw
RNApolII=/Volumes/external.data/MeisterLab/bisiaka/supercoilingdata/supercoiling2/fountainspaperfiles/bigwigs/RNAPII_TIR1_noauxin_mean.bw
TOP1=/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP1-GFP_AM05_L2-L3_ce11_avr.bw
TOP2=/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP2-GFP_MDX53_L2-L3_ce11_avr.bw

computeMatrix scale-regions -S $COH1 $RNApolII $TOP1 $TOP2 \
  -R fountains_assymmetryQ1.bed fountains_assymmetryQ3.bed fountains_assymmetryQ5.bed \
  -b 20000 -a 20000  --regionBodyLength 2000 \
  -o matrix_asymQ_COH1all.gz  --binSize 100 \
  --outFileNameMatrix matrix_asymQ_COH1all.tab \
  --outFileSortedRegions matrix_asymQ_COH1all.bed \
  --startLabel "fountain" --endLabel "tip" \
  --sortRegions descend  --samplesLabel "COH-1" "RNApolII" "TOP-1" "TOP-2" \
  -p 2


plotHeatmap -m matrix_asymQ_COH1all.gz \
      -out matrix_asymQ_COH1_heatmap_all.png \
      --regionsLabel "AsymQ1" "AsymQ3" "AsymQ5" \
      --plotTitle "Asymmetric quantiles" \
      --colorMap Reds Oranges Blues Purples Greens Greys\
      --startLabel "fountain" --endLabel "tip" \
      --zMin 0 -1.5 -1.5 -0.4  \
      --zMax 2 1.5 1.5 2  \
      --yMin 0.5 -1 -1 -1 \
      --yMax 2 5 3 3.5 \
      --linesAtTickMarks






plotProfile -m matrix_asymQ_COH1all.gz \
              -out matrix_asymQ_COH1_profile_all.png \
              --numPlotsPerRow 4 \
              --regionsLabel "AsymQ1" "AsymQ3" "AsymQ5" \
              --startLabel "fountain" --endLabel "tip" \
              --yMin 0.5 -1 -1 -1 \
              --yMax 2 5 3 3.5 \
              --colors "cyan" "black" "magenta" \
              --plotTitle "Asymmetric quantiles" 
              