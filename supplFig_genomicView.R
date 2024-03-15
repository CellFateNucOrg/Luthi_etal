library(plotgardener)
library(BSgenome.Celegans.UCSC.ce11)
library(HiCExperiment)
library(GenomicRanges)
library(HiContacts)
library(RColorBrewer)
library(rtracklayer)
library(dplyr)
library(TxDb.Celegans.UCSC.ce11.refGene)

projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")


gr<-GRanges(seqnames="chrII",ranges=IRanges(start=9181001, end=9381000))
seqlevels(gr)<-seqlevels(Celegans)

# Data -----
hic828<-"/Volumes/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/cohesin1kb_files/mcool/828_merge_1000.mcool"
cf828<-CoolFile(hic828)
cool828 <- import(cf828,focus="chrII:9000000-9500000",resolution=1000)

hic366<-"/Volumes/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/cohesin1kb_files/mcool/366_merge_1000.mcool"
cf366<-CoolFile(hic366)
cool366 <- import(cf366,focus="chrII:9000000-9500000",resolution=1000)


div_contacts <- divide(cool828, by = cool366)
div_contacts1<-subsetByOverlaps(div_contacts,gr)

hicplotted<-plotMatrix(
  div_contacts1,
  use.scores = 'balanced.fc',
  scale = 'log2',
  limits = c(-1, 1),
  cmap = bwrColors(), maxDistance=70000,
  caption=F
)

hicplotted

# plotGG(
#   plot = hicplotted,
#   x = 0, y = 0,
#   width = 18, height = 4, default.units = "cm",
#   just = c("left", "top")
# )

## moving to plotgardener framework
mat<-data.frame(chrII=start(anchors(div_contacts1)$first),chrII=start(anchors(div_contacts1)$second),count=scores(div_contacts1)$balanced.l2fc)


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
seqlevels(fountains)<-seqlevels(Celegans)

daughertyL3<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))
daughertyL3<-daughertyL3[daughertyL3$L3_chromHMMState %in% c("L3_activeEnhancer")]

JaenesL3<-readRDS(paste0(publicDataDir,"/Jaenes2018_enhancers_ce11_stages_chromHMM.rds"))
JaenesL3<-JaenesL3[JaenesL3$topState_L3_chromHMM %in% c("Active enhancer")]

daughertyATAC<-import("/Volumes/external.data/MeisterLab/publicData/ATACseq/Daugherty2017_GenomeRes_GSE89608/atac_PRJNA352701/ATAC_L3_avr_Daugherty2017.bw",
                      selection=grFix)

jaenesATAC<-import("/Volumes/external.data/MeisterLab/publicData/ATACseq/Jaenes2018_eLife_GSE114494/GSE114439_atac_wt_l3_ce11.bw",
                   selection=grFix)

coh1<-import("/Volumes/MeisterLab/publicData/ChIPseq_realigned/modEncode_SMC/chip-seq-pipeline2/bigwig_fc/COH1_fem2_AD_SDQ0809_rep.pooled_x_ctl.pooled.fc.signal.bigwig",
             selection=grFix)

H3K27ac<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/GSM624432_WA30634849_H3K27AC_N2_L3_1_ce11.bw",
                selection=grFix)

#H3K27me3<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/GSM562734_HK00013_H3K27ME31E7_N2_L3_1_ce11.bw",
#                 selection=grFix)

H3K27me3<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_realigned/modEncode_histone/bigwig_fc/H3K27me3_UP07449_rep.pooled_x_ctl.pooled.fc.signal.bigwig",
                 selection=grFix)


top1<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP1-GFP_AM05_L2-L3_ce11_avr.bw",
             selection=grFix)

top2<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP2-GFP_MDX53_L2-L3_ce11_avr.bw",
             selection=grFix)

chromStates<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromStates_L3_Evans2016_ce11_rgb.bed")

chromDomains<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromDomains_L3_Evans2016_ce11_rgb.bed")

chromHMM<-import("/Volumes/external.data/MeisterLab/publicData/ATACseq/Daugherty2017_GenomeRes_GSE89608/chromHMM/chromHMM_L3_Daugherty2017_ce11_rgb.bed")

tpm366<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/tracks/PMW366_TPM_avr.bw",
               selection=grFix)

coh1lfc<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/tracks/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_lfc.bw",
                selection=grFix)



#change gr to fit hic matrix
grFix<-GRanges(seqnames=seqnames(gr),ranges=IRanges(start=min(mat$chrII),end=max(mat$chrII.1)))
seqlevels(grFix)<-seqlevels(Celegans)




############ plots -----
params <- pgParams(
  chrom = as.character(seqnames(grFix)), chromstart = start(grFix), chromend = end(grFix),
  assembly = "ce11",
  just = c("left", "top"),
  default.units = "cm",
  x=0.5,
  width =17,
  fontsize = 10,
  spaceHeight=0, spaceWidth=0
)

pageCreate(width = 18, height = 20,params=params)
hicHeight=5
chipBin=200
hicPlot<-plotHicTriangle(
  data = mat,
  params=params,
  palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  zrange=c(-1,1),
  x = 0.5, y = 0, width =17, height=hicHeight,
  colorTrans="linear"
)

annoHeatmapLegend(
  plot = hicPlot,
  params=params,
  x = 17, y = 0.2, width = 0.3, height = 2,
  fontcolor="black"
)

plotText(label = "HiC\nlog2(COH-1/TEVonly)", fontcolor = "black",
          y = 0.2, params=params)

annoGenomeLabel(
  plot = hicPlot,
  params=params,
  y = hicHeight, height = 0.25,
)

## fountains -----
p<-plotRanges(fountains,
           params=params,
           y = hicHeight+0.6, height=0.3,
           boxHeight=0.3,
           fill="red")
plotText(label = "Fountain tip", fontcolor = "red",
  y = hicHeight+0.6, params=params)


## enahncers -----
p<-plotRanges(daughertyL3,
           params=params,
           y = hicHeight+1.1, height=0.3,
           boxHeight=0.3,
           fill="violet")
plotText(label = "Enhancers(Daugherty)", fontcolor = "violet",
  x = 0.5, y = hicHeight+1.1, params=params)

p<-plotRanges(JaenesL3,
           params=params,
           y = hicHeight+1.6, height=0.3,
           spaceHeight=0,spaceWidth=0,
           boxHeight=0.3,
           fill="violet")
plotText(label = "Enhancers(Jaenes)", fontcolor = "violet",
         y = hicHeight+1.6, params=params)

## ATAC and chip seq----
p<-plotSignal(coh1,
           params=params,
           binSize=chipBin,
           fill="red2",
           y = hicHeight+2, height=1,
           linecolor="red2")
plotText(label = "COH-1 ChIP", fontcolor = "red2",
         y = hicHeight+2.1, params=params)


p<-plotSignal(daughertyATAC,
           params=params,
           binSize=chipBin,
           fill="violet",
           y = hicHeight+3, height=1,
           linecolor="violet")
plotText(label = "ATAC(Daugherty)", fontcolor = "violet",
         y = hicHeight+3.1, params=params)


p<-plotSignal(jaenesATAC,
           params=params,
           binSize=chipBin,
           fill="violet",
           y = hicHeight+4, height=1,
           linecolor="violet")
plotText(label = "ATAC(Jaenes)", fontcolor = "violet",
         y = hicHeight+4.1, params=params)

p<-plotSignal(H3K27ac,
           params=params,
           binSize=chipBin,
           fill="darkgreen",
           y = hicHeight+5, height=1,
           linecolor="darkgreen",
           range=c(0,max(H3K27ac$score)),negData=T)
plotText(label = "H3K27Ac ChIP", fontcolor = "darkgreen",
         y = hicHeight+5.1, params=params)

p<-plotSignal(H3K27me3,
           params=params,
           binSize=chipBin,
           fill="darkgrey",
           y = hicHeight+6, height=1,
           linecolor="darkgrey",
           range=c(0,max(H3K27me3$score)))
plotText(label = "H3K27me3 ChIP", fontcolor = "darkgrey",
         y = hicHeight+6.1, params=params)

p<-plotSignal(top1,
           params=params,
           binSize=chipBin,
           fill="orange2",
           y = hicHeight+7, height=1,
           linecolor="orange2",
           negData=T)
plotText(label = "TOP-1 ChIP", fontcolor = "orange2",
         y = hicHeight+7.1, params=params)

p<-plotSignal(top2,
           params=params,
           binSize=chipBin,
           fill="orange4",
           y = hicHeight+8, height=1,
           linecolor="orange4",
           negData=T)
plotText(label = "TOP-2 ChIP", fontcolor = "orange4",
         y = hicHeight+8.1, params=params)

## chromHMM -----
p<-plotRanges(chromDomains,
           params=params,
           y = hicHeight+9.4, height=0.3,
           boxHeight=0.3, fill=chromDomains$itemRgb)
plotText(label = "ChromatinDomains (Evans) ", fontcolor = "black",
         y = hicHeight+9.1, params=params)


p<-plotRanges(chromStates,
           params=params,
           y = hicHeight+10.1, height=0.3,
           boxHeight=0.3,fill=chromStates$itemRgb)
plotText(label = "ChromatinStates (Evans)", fontcolor = "black",
         y = hicHeight+9.8, params=params)


p<-plotRanges(chromHMM,
           params=params,
           y = hicHeight+10.8, height=0.3,
           boxHeight=0.3,fill=chromHMM$itemRgb)
plotText(label = "ChromainHMM (Daugherty)", fontcolor = "black",
         y = hicHeight+10.5, params=params)

### RNAseq -----
p<-plotSignal(coh1lfc,
           params=params,
           binSize=chipBin,
           fill="red2",
           y = hicHeight+11.2, height=1,
           linecolor="red2",
           negData=T)
plotText(label = "COH-1 cleavage RNA log2FC", fontcolor = "red2",
         y = hicHeight+11.1, params=params)

p<-plotSignal(tpm366,
           params=params,
           binSize=chipBin,
           fill="blue4",
           y = hicHeight+12.2, height=1,
           linecolor="blue4",
           negData=T)
plotText(label = "TEV only TPM", fontcolor = "blue4",
         y = hicHeight+12.3, params=params)


p<-plotGenes(params=params,
                y=hicHeight+13.3, height=1,
                fontsize=1)

fountLoc<-subsetByOverlaps(fountains,grFix)

annoHighlight(
  plot=p,
  params=params,
  chrom = as.character(seqnames(fountLoc)),
  chromstart = start(fountLoc), chromend = end(fountLoc),
  y = hicHeight, height = 14.3
)

pageGuideHide()

