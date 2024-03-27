library(plotgardener)
library(BSgenome.Celegans.UCSC.ce11)
library(HiCExperiment)
library(GenomicRanges)
library(HiContacts)
library(RColorBrewer)
library(rtracklayer)
library(dplyr)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(InteractionSet)

projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
fountainFigDir=paste0(projectDir,"/fountainFigures")

fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
seqlevels(fountains)<-seqlevels(Celegans)

resolution=2000
gr<-GRanges(seqnames="chrIV",ranges=IRanges(start=5550000, end=5860000))
seqlevels(gr)<-seqlevels(Celegans)

# Data -----
hic828<-"/Volumes/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/cohesin1kb_files/mcool/828_merge_1000.mcool"
cf828<-CoolFile(hic828)
cool828 <- import(cf828,focus=paste0(seqnames(gr),":",start(gr),"-",end(gr)),
                  resolution=resolution)

hic366<-"/Volumes/external.data/MeisterLab/mdas/Illumina_HiC/2021_HiCworked/cohesin1kb_files/mcool/366_merge_1000.mcool"
cf366<-CoolFile(hic366)
cool366 <- import(cf366,focus=paste0(seqnames(gr),":",start(gr),"-",end(gr)),
                  resolution=resolution)


div_contacts <- divide(cool828, by = cool366)
div_contacts1<-subsetByOverlaps(div_contacts,gr)


## moving to plotgardener framework
mat366<-data.frame(chr=start(anchors(cool366)$first),chr=start(anchors(cool366)$second),count=scores(cool366)$balanced)
mat828<-data.frame(chr=start(anchors(cool828)$first),chr=start(anchors(cool828)$second),count=scores(cool828)$balanced)
mat<-data.frame(chr=start(anchors(div_contacts1)$first),chr=start(anchors(div_contacts1)$second),count=scores(div_contacts1)$balanced.l2fc)


#change gr to fit hic matrix
#grFix<-GRanges(seqnames=seqnames(gr),ranges=IRanges(start=min(mat$chr),end=max(mat$chr.1)))
grFix<-gr
seqlevels(grFix)<-seqlevels(Celegans)

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


chromDomains<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromDomains_L3_Evans2016_ce11_rgb.bed")



############ plots -----

pdf(file=paste0(finalFigDir,"/fig2ef_skn-1_square.pdf"),width=5,height=16.5,paper="a4")

params <- pgParams(
  chrom = as.character(seqnames(grFix)), chromstart = start(grFix), chromend = end(grFix),
  assembly = "ce11",
  just = c("left", "top"),
  default.units = "cm",
  x=0.5,
  width =4.5,
  fontsize = 10,
  spaceHeight=0, spaceWidth=0
)

textHeight=0.2
boxHeight=0.25
rangesHeight=0.25
signalHeight=1
hicHeight=4.5
transcriptHeight=3
chipBin=200
geneHeight=0.8

currentHeight=0
pageCreate(width = 5.5, height = 16.5, params=params, showGuides=FALSE)

hicPlot<-plotHicSquare(
  data = mat366,
  params=params,
  palette=colorRampPalette(colors=c("white","#4a90c2","black","white"),interpolate="spline"),
  zrange=c(1e-4,1e-2),
  x = 0.5, y = currentHeight, height=hicHeight,
  colorTrans="log10"
)

annoHeatmapLegend(
  plot = hicPlot,
  params=params,
  x = 0, y = currentHeight+hicHeight/2, width = 0.3, height = 2,
  fontcolor="black",
  just=c("left","center"),
  scientific=T,font=7
)

plotText(label = "TEV control", fontcolor = "white",
         x=0.6, y = currentHeight+0.1, params=params,fontface="bold")
currentHeight<-currentHeight+hicHeight+0.1


hicPlot<-plotHicSquare(
  data = mat828,
  params=params,
  palette=colorRampPalette(colors=c("white","#4a90c2","black","white"),interpolate="spline"),
  zrange=c(1e-4,1e-2),
  x = 0.5, y = currentHeight, height=hicHeight,
  colorTrans="log10"
)

annoHeatmapLegend(
  plot = hicPlot,
  params=params,
  x = 0, y = currentHeight+hicHeight/2, width = 0.3, height = 2,
  fontcolor="black",
  just=c("left","center"),
  scientific=T,font=7
)

plotText(label =bquote("cohesin"^"COH-1"~"cleavage"), fontcolor = "white",
         x=0.6,y = currentHeight+0.1, params=params, fontface="bold")
currentHeight<-currentHeight+hicHeight+0.2


hicPlot<-plotHicSquare(
  data = mat,
  params=params,
  palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  zrange=c(-1,1),
  x = 0.5, y = currentHeight, height=hicHeight,
  colorTrans="linear"
)

annoHeatmapLegend(
  plot = hicPlot,
  params=params,
  x = 0, y = currentHeight+hicHeight/2, width = 0.3, height = 2,
  fontcolor="black",
  just=c("left","center"),
  scientific=F,font=7
)

plotText(label = bquote("cohesin"^"COH-1"~"/TEV control"), fontcolor = "black",
          x=0.6,y = currentHeight+0.1, params=params,fontface="bold")

currentHeight<-currentHeight+ hicHeight

## fountains -----
p<-plotRanges(fountains,
           params=params,
           y = currentHeight, height=rangesHeight,
           boxHeight=boxHeight,
           fill="red")

currentHeight<-currentHeight+rangesHeight + 0.05


## enahncers -----
p<-plotRanges(daughertyL3,
           params=params,
           y = currentHeight, height=rangesHeight,
           boxHeight=boxHeight,
           spaceHeight=0,spaceWidth=0,
           fill="violet")
currentHeight<-currentHeight+rangesHeight + 0.05

p<-plotRanges(JaenesL3,
           params=params,
           y = currentHeight, height=rangesHeight,
           spaceHeight=0,spaceWidth=0,
           boxHeight=boxHeight,
           fill="violet")
currentHeight<-currentHeight+rangesHeight + 0.05

p<-plotGenes(params=params,
             y=currentHeight-0.2, height=signalHeight,
             fontsize=0,
             geneHighlights=data.frame(gene="skn-1",color="black"))
currentHeight<-currentHeight + signalHeight/2
plotText(label = "skn-1", fontcolor = "black",
         y = currentHeight, params=params,fontface="italic",
         x=1.5)
currentHeight<-currentHeight+textHeight+0.1


# p<-plotSignal(coh1,
#               params=params,
#               binSize=chipBin,
#               fill="red2",
#               y = currentHeight, height=signalHeight,
#               linecolor="red2")
# currentHeight<-currentHeight + signalHeight + 0.1


## chromHMM -----
p<-plotRanges(chromDomains,
              params=params,
              y = currentHeight, height=rangesHeight,
              boxHeight=boxHeight, fill=chromDomains$itemRgb)
currentHeight<-currentHeight +rangesHeight


annoGenomeLabel(
  plot = hicPlot,
  params=params,
  y = currentHeight, height = textHeight,
  scale="Mb",fontsize=8
)

dev.off()


