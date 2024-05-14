#! /bin/bash
micromamba activate deeptools

COH1=/Volumes/external.data/MeisterLab/bisiaka/supercoilingdata/supercoiling2/fountainspaperfiles/bigwigs/COH-1modencode_mean.bw
RNApolII=/Volumes/external.data/MeisterLab/bisiaka/supercoilingdata/supercoiling2/fountainspaperfiles/bigwigs/RNAPII_TIR1_noauxin_mean.bw
TOP1=/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP1-GFP_AM05_L2-L3_ce11_avr.bw
TOP2=/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP2-GFP_MDX53_L2-L3_ce11_avr.bw
H3K27ac=/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/GSM624432_WA30634849_H3K27AC_N2_L3_1_ce11.bw
H3K27me3=/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/GSM562734_HK00013_H3K27ME31E7_N2_L3_1_ce11.bw


computeMatrix scale-regions -S $COH1 $RNApolII $TOP1 $TOP2 $H3K27ac $H3K27me3 \
  -R fountains_assymmetryQ1.bed fountains_assymmetryQ3.bed fountains_assymmetryQ5.bed \
  -b 20000 -a 20000  --regionBodyLength 2000 \
  -o matrix_asymQ_COH1all.gz  --binSize 100 \
  --outFileNameMatrix matrix_asymQ_COH1all.tab \
  --outFileSortedRegions matrix_asymQ_COH1all.bed \
  --startLabel "fountain" --endLabel "tip" \
  --sortRegions descend  --samplesLabel "COH-1" "RNApolII" "TOP-1" "TOP-2" "H3K27ac" "H3K27me3" \
  -p 2


plotHeatmap -m matrix_asymQ_COH1all.gz \
      -out matrix_asymQ_COH1_heatmap_all.png \
      --regionsLabel "AsymQ1" "AsymQ3" "AsymQ5" \
      --plotTitle "Asymmetric quantiles" \
      --colorMap Reds Oranges Blues Purples Greens Greys\
      --startLabel "fountain" --endLabel "tip" \
      --zMin 0 -1.5 -1.5 -0.4 -1 -1\
      --zMax 2 1.5 1.5 2  2 1\
      --yMin 0.5 -1 -1 -1 -0.5 -1 \
      --yMax 2 5 3 3.5 1.5 0.5 \
      --linesAtTickMarks






plotProfile -m matrix_asymQ_COH1all.gz \
              -out matrix_asymQ_COH1_profile_all.png \
              --numPlotsPerRow 6 \
              --regionsLabel "AsymQ1" "AsymQ3" "AsymQ5" \
              --startLabel "fountain" --endLabel "tip" \
              --yMin 0.5 -1 -1 -1 -0.5 -1 \
              --yMax 2 5 3 3.5 1.5 0.5 \
              --colors "cyan" "black" "magenta" \
              --plotTitle "Asymmetric quantiles"
