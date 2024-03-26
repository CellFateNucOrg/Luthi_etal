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
i=500
resolution=2000
#for(i in 1:length(fountains)){
gr<-GRanges(seqnames="chrIV",ranges=IRanges(start=5584000, end=6110000))
skn1zoom<-GRanges(seqnames="chrIV",ranges=IRanges(start=5644000, end=5676000))
nhr46zoom<-GRanges(seqnames="chrIV",ranges=IRanges(start=5744000, end=5768000))

#seqlevels(gr)<-seqlevels(Celegans)
#gr<-resize(fountains[i],width=max(regions)+resolution,fix="center")

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
mat<-data.frame(chr=start(anchors(div_contacts1)$first),chr=start(anchors(div_contacts1)$second),count=scores(div_contacts1)$balanced.l2fc)


#change gr to fit hic matrix
grFix<-GRanges(seqnames=seqnames(gr),ranges=IRanges(start=min(mat$chr),end=max(mat$chr.1)))
seqlevels(grFix)<-seqlevels(Celegans)

daughertyL3<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))
daughertyL3<-daughertyL3[daughertyL3$L3_chromHMMState %in% c("L3_activeEnhancer")]

JaenesL3<-readRDS(paste0(publicDataDir,"/Jaenes2018_enhancers_ce11_stages_chromHMM.rds"))
JaenesL3<-JaenesL3[JaenesL3$topState_L3_chromHMM %in% c("Active enhancer")]

daughertyATAC<-import("/Volumes/external.data/MeisterLab/publicData/ATACseq/Daugherty2017_GenomeRes_GSE89608/atac_PRJNA352701/ATAC_L3_avr_Daugherty2017.bw",
                      selection=grFix)

jaenesATAC<-import("/Volumes/external.data/MeisterLab/publicData/ATACseq/Jaenes2018_eLife_GSE114494/GSE114439_atac_wt_l3_ce11.bw",
                   selection=grFix)

intestineATAC <- import("/Volumes/MeisterLab/publicData/ATACseq/Durham2021_L2tissueATAC_GSE157017/intestine.bigWig", selection=grFix)

neuronATAC <- import("/Volumes/MeisterLab/publicData/ATACseq/Durham2021_L2tissueATAC_GSE157017/neuron.bigWig", selection=grFix)

coh1<-import("/Volumes/MeisterLab/publicData/ChIPseq_realigned/modEncode_SMC/chip-seq-pipeline2/bigwig_fc/COH1_fem2_AD_SDQ0809_rep.pooled_x_ctl.pooled.fc.signal.bigwig",
             selection=grFix)


chromStates<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromStates_L3_Evans2016_ce11_rgb.bed")

chromDomains<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromDomains_L3_Evans2016_ce11_rgb.bed")


deletions <- import(paste0(fountainsDir,"/ubs_deletion.bed"))
deletions$name<-gsub("_deletion","",deletions$name)

arcC <- import("/Volumes/MeisterLab/publicData/InteractionData/Huang-Ahringer_ARC-C_GSE144673/Huang_skn1SignificantInteractions.bedpe")
arcCgi<-makeGInteractionsFromGRangesPairs(arcC)

############ plots -----

pdf(file=paste0(finalFigDir,"/supplFig_genomicView_skn1deletions.pdf"),width=18,height=29,paper="a4")

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

textHeight=0.2
boxHeight=0.3
rangesHeight=0.4
signalHeight=1
hicHeight=4
transcriptHeight=2.5
chipBin=200

currentHeight=0
pageCreate(width = 18, height = 29, params=params, showGuides=FALSE)

plotText(label = "Isiaka et al., Supl. Figure", fontcolor = "black",
         y = currentHeight+1, params=params, x=9, font=14)

currentHeight<-currentHeight+2

hicPlot<-plotHicTriangle(
  data = mat,
  params=params,
  palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  zrange=c(-1,1),
  y = currentHeight, width =17, height=hicHeight,
  colorTrans="linear"
)

annoHeatmapLegend(
  plot = hicPlot,
  params=params,
  x = 17, y = currentHeight+0.3, width = 0.3, height = 2,
  fontcolor="black"
)

plotText(label = "HiC", fontcolor = "black",
          y = currentHeight+0.3, params=params)

currentHeight<-currentHeight+ hicHeight

hicPlot<-plotPairsArches(
  data = arcCgi, params = params,
  fill = colorby("score", palette =
                   colorRampPalette(c("dodgerblue3", "darkblue"))),
  linecolor = "fill",
  archHeight = "score",
  alpha = 0.7,
  y = currentHeight, height = 1.2,
  flip=TRUE
)
plotText(label = "ArcC", fontcolor = "darkblue",
         y = currentHeight+1.2, params=params)
currentHeight<-currentHeight+ 1.5 + textHeight + 0.1

## fountains -----
hicPlot<-plotRanges(fountains,
           params=params,
           y = currentHeight+0.2, height=rangesHeight,
           boxHeight=boxHeight,
           fill="red")
plotText(label = "Fountain tip", fontcolor = "red",
         y = currentHeight, params=params)
currentHeight<-currentHeight+rangesHeight+textHeight+0.1

## enahncers -----
hicPlot<-plotRanges(daughertyL3,
                     params=params,
                     y = currentHeight+textHeight, height=rangesHeight,
                     boxHeight=boxHeight,
                     spaceHeight=0,spaceWidth=0,
                     fill="violet")
plotText(label = "Enhancers(Daugherty)", fontcolor = "violet",
                   y = currentHeight, params=params)
currentHeight<-currentHeight+rangesHeight +textHeight +0.1

hicPlot<-plotRanges(JaenesL3,
                     params=params,
                     y = currentHeight+textHeight, height=rangesHeight,
                     spaceHeight=0,spaceWidth=0,
                     boxHeight=boxHeight,
                     fill="violet")

plotText(label = "Enhancers(Jaenes)", fontcolor = "violet",
                   y = currentHeight, params=skn1params)
currentHeight<-currentHeight+rangesHeight+ textHeight+ 0.1


hicPlot<-plotSignal(coh1,
              params=params,
              binSize=chipBin,
              fill="red2",
              y = currentHeight, height=signalHeight,
              linecolor="red2")
plotText(label = "COH-1 ChIP", fontcolor = "red2",
         y = currentHeight, params=params)
currentHeight<-currentHeight + signalHeight + 0.1


## chromHMM -----
hicPlot<-plotText(label = "ChromatinDomains (Evans) ", fontcolor = "black",
         y = currentHeight, params=params)
hicPlot<-plotRanges(chromDomains,
              params=params,
              y = currentHeight+textHeight, height=rangesHeight,
              boxHeight=boxHeight, fill=chromDomains$itemRgb)
currentHeight<-currentHeight + textHeight+rangesHeight

hicPlot<-plotGenes(params=params,
             y=currentHeight, height=signalHeight,
             fontsize=0, fill = c("#669fd9", "#abcc8e"),
             geneHighlights=data.frame(gene="skn-1",color="darkgreen"))
plotText(label = "skn-1", fontcolor = "darkgreen",
                  y = currentHeight+signalHeight-0.3, params=params,
                  x=2.6)

currentHeight<-currentHeight + signalHeight + textHeight


# fountLoc<-subsetByOverlaps(fountains,grFix)
#
# for(j in 1:length(fountLoc)){
#   annoHighlight(
#     plot=hicPlot,
#     params=params,
#     chrom = as.character(seqnames(fountLoc[j])),
#     chromstart = start(fountLoc[j]), chromend = end(fountLoc[j]),
#     y = hicHeight, height = currentHeight-hicHeight
#   )
# }
annoHighlight(
  plot=hicPlot,
  params=params,
  chrom = as.character(seqnames(skn1zoom)),
  chromstart = start(skn1zoom), chromend = end(skn1zoom),
  y = hicHeight+2, height = currentHeight-hicHeight-2
)

annoHighlight(
  plot=hicPlot,
  params=params,
  chrom = as.character(seqnames(nhr46zoom)),
  chromstart = start(nhr46zoom), chromend = end(nhr46zoom),
  y = hicHeight+2, height = currentHeight-hicHeight-2
)
#pageGuideHide()

annoGenomeLabel(
  plot = hicPlot,
  params=params,
  y = currentHeight, height = textHeight,
  margin=-0.3
)


skn1params <- pgParams(
  chrom = as.character(seqnames(skn1zoom)), chromstart = start(skn1zoom), chromend = end(skn1zoom),
  assembly = "ce11",
  just = c("left", "top"),
  default.units = "cm",
  x=0.5,
  width =10,
  fontsize = 10,
  spaceHeight=0, spaceWidth=0
)

annoZoomLines(
  plot = hicPlot, params=skn1params,linecolor="black",
  y0 = currentHeight, x1 = c(0.5, 10.5), y1 = currentHeight+0.5,
  lty=1
)

nhr46params <- pgParams(
  chrom = as.character(seqnames(nhr46zoom)), chromstart = start(nhr46zoom), chromend = end(nhr46zoom),
  assembly = "ce11",
  just = c("left", "top"),
  default.units = "cm",
  x=11,
  width =6.5,
  fontsize = 10,
  spaceHeight=0, spaceWidth=0
)

annoZoomLines(
  plot = hicPlot, params=nhr46params,linecolor="black",
  y0 = currentHeight, x1 = c(11, 17.5), y1 = currentHeight+0.5,
  lty=1
)
currentHeight<-currentHeight + 0.5

######## subplots -----


skn1Plot<-plotHicTriangle(
  data = mat,
  params=skn1params,
  palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  zrange=c(-1,1),
  y = currentHeight, height=hicHeight-2,
  colorTrans="linear"
)
nhr46Plot<-plotHicTriangle(
  data = mat,
  params=nhr46params,
  palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  zrange=c(-1,1),
  y = currentHeight, height=hicHeight-2,
  colorTrans="linear"
)
plotText(label = "HiC", fontcolor = "black",
         y = currentHeight+0.2, params=skn1params)
currentHeight<-currentHeight + hicHeight-2

subplotHeight<-currentHeight

skn1Plot<-plotPairsArches(
  data = arcCgi, params = skn1params,
  fill = colorby("score", palette =
                   colorRampPalette(c("dodgerblue3", "darkblue"))),
  linecolor = "fill",
  archHeight = "score",
  alpha = 0.4, clip=F, clip.noanchor=F,
  y = currentHeight, height = 1.2, linewidth=1,
  flip=TRUE
)

nhr46Plot<-plotPairsArches(
  data = arcCgi, params = nhr46params,
  fill = colorby("score", palette =
                   colorRampPalette(c("dodgerblue3", "darkblue"))),
  linecolor = "fill",
  archHeight = "score",
  alpha = 0.7, clip=F, clip.noanchor=F,
  y = currentHeight, height = 1.2,
  flip=TRUE, curvature=15
)

skn1Plot<-plotText(label = "ArcC", fontcolor = "darkblue",
                  y = currentHeight, params=skn1params)
currentHeight<-currentHeight+ 1.2


skn1Plot<-plotRanges(chromStates,
              params=skn1params,
              y = currentHeight+textHeight, height=rangesHeight,
              boxHeight=boxHeight,fill=chromStates$itemRgb)
nhr46Plot<-plotRanges(chromStates,
                     params=nhr46params,
                     y = currentHeight+textHeight, height=rangesHeight,
                     boxHeight=boxHeight,fill=chromStates$itemRgb)
skn1Plot<-plotText(label = "ChromatinStates (Evans)", fontcolor = "black",
                   y = currentHeight, params=skn1params)
currentHeight<-currentHeight+rangesHeight + textHeight


## fountains -----
skn1Plot<-plotRanges(fountains,
                    params=skn1params,
                    y = currentHeight, height=rangesHeight,
                    boxHeight=boxHeight,
                    fill="red")
nhr46Plot<-plotRanges(fountains,
                    params=nhr46params,
                    y = currentHeight, height=rangesHeight,
                    boxHeight=boxHeight,
                    fill="red")
skn1Plot<-plotText(label = "Fountain tip", fontcolor = "red",
                  y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight+rangesHeight


## enahncers -----
skn1Plot<-plotRanges(daughertyL3,
                    params=skn1params,
                    y = currentHeight, height=rangesHeight,
                    boxHeight=boxHeight,
                    spaceHeight=0,spaceWidth=0,
                    fill="violet")
nhr46Plot<-plotRanges(daughertyL3,
                     params=nhr46params,
                     y = currentHeight, height=rangesHeight,
                     boxHeight=boxHeight,
                     spaceHeight=0,spaceWidth=0,
                     fill="violet")
skn1Plot<-plotText(label = "Enhancers(Daugherty)", fontcolor = "violet",
                  y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight+rangesHeight

skn1Plot<-plotRanges(JaenesL3,
                    params=skn1params,
                    y = currentHeight, height=rangesHeight,
                    spaceHeight=0,spaceWidth=0,
                    boxHeight=boxHeight,
                    fill="violet")
nhr46Plot<-plotRanges(JaenesL3,
                     params=nhr46params,
                     y = currentHeight, height=rangesHeight,
                     spaceHeight=0,spaceWidth=0,
                     boxHeight=boxHeight,
                     fill="violet")
skn1Plot<-plotText(label = "Enhancers(Jaenes)", fontcolor = "violet",
                  y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight+rangesHeight

# ATAC L3
skn1Plot<-plotSignal(daughertyATAC,
                    params=skn1params,
                    binSize=chipBin/10,
                    fill="violet",
                    y = currentHeight, height=signalHeight,
                    linecolor="violet")
nhr46Plot<-plotSignal(daughertyATAC,
                     params=nhr46params,
                     binSize=chipBin/10,
                     fill="violet",
                     y = currentHeight, height=signalHeight,
                     linecolor="violet")
plotText(label = "ATAC(Daugherty)", fontcolor = "violet",
         y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight + signalHeight


skn1Plot<-plotSignal(jaenesATAC,
                    params=skn1params,
                    binSize=chipBin/10,
                    fill="violet",
                    y = currentHeight, height=signalHeight,
                    linecolor="violet")
nhr46Plot<-plotSignal(jaenesATAC,
                    params=nhr46params,
                    binSize=chipBin/10,
                    fill="violet",
                    y = currentHeight, height=signalHeight,
                    linecolor="violet")
plotText(label = "ATAC(Jaenes)", fontcolor = "violet",
         y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight + signalHeight



## ATAC tissue----
skn1Plot<-plotSignal(intestineATAC,
                     params=skn1params,
                     binSize=chipBin/10,
                     fill="yellow3",
                     y = currentHeight, height=signalHeight,
                     linecolor="yellow3")
nhr46Plot<-plotSignal(intestineATAC,
                    params=nhr46params,
                    binSize=chipBin/10,
                    fill="yellow3",
                    y = currentHeight, height=signalHeight,
                    linecolor="yellow3")
plotText(label = "Intestine ATAC", fontcolor = "yellow3",
                  y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight + signalHeight


skn1Plot<-plotSignal(neuronATAC,
                    params=skn1params,
                    binSize=chipBin/10,
                    fill="purple3",
                    y = currentHeight, height=signalHeight,
                    linecolor="purple3")
nhr46Plot<-plotSignal(neuronATAC,
                     params=nhr46params,
                     binSize=chipBin/10,
                     fill="purple3",
                     y = currentHeight, height=signalHeight,
                     linecolor="purple3")
plotText(label = "Neuron ATAC", fontcolor = "purple3",
                  y = currentHeight+0.1, params=skn1params)
currentHeight<-currentHeight + signalHeight

skn1Plot<-plotRanges(deletions,
              params=skn1params,
              y = currentHeight, height=rangesHeight*3,
              spaceHeight=0.1,spaceWidth=0,
              boxHeight=boxHeight,
              fill="black")
nhr46Plot<-plotRanges(deletions,
                     params=nhr46params,
                     y = currentHeight, height=rangesHeight*3,
                     spaceHeight=0.1,spaceWidth=0,
                     boxHeight=boxHeight,
                     fill="black")
plotText(label = "deletions", fontcolor = "black",
         y = currentHeight+0.1, params=skn1params)
plotText(label=sort(deletions)$name[1],y = currentHeight+rangesHeight*2+0.1,x=2.6,
         fontcolor="black",params=skn1params)
plotText(label=sort(deletions)$name[2],y = currentHeight+rangesHeight*2+0.1,x=5.5,
         fontcolor="black",params=skn1params)
plotText(label=sort(deletions)$name[3],y = currentHeight+rangesHeight*1+0.1,x=5.5,
         fontcolor="black",params=skn1params)
plotText(label=sort(deletions)$name[4],y = currentHeight+0.2,x=5,
         fontcolor="black",params=skn1params)
plotText(label=sort(deletions)$name[5],y = currentHeight+rangesHeight*2+0.1,x=8,
         fontcolor="black",params=skn1params)
plotText(label=sort(deletions)$name[6],y = currentHeight+rangesHeight*1+0.1,x=8,
         fontcolor="black",params=skn1params)

plotText(label=sort(deletions)$name[7],y = currentHeight+rangesHeight*1+0.1,x=13.5,
         fontcolor="black",params=nhr46params)

currentHeight<-currentHeight+rangesHeight*3+0.1

skn1Plot<-plotTranscripts(params=skn1params,
                   y=currentHeight, height=transcriptHeight, labels="gene",
                   fontsize=9,transcriptHighlights=data.frame(gene="skn-1",color="darkgreen"))
nhr46Plot<-plotTranscripts(params=nhr46params,
                          y=currentHeight, height=transcriptHeight, labels="gene",
                          fontsize=9)
currentHeight<-currentHeight + transcriptHeight + 0.1

reducedDel<-reduce(deletions)
for(j in 1:length(reducedDel)){
  if(start(reducedDel[j])<end(skn1zoom)){
    annoHighlight(
      plot=skn1Plot,
      params=skn1params,
      chrom = as.character(seqnames(reducedDel[j])),
      chromstart = start(reducedDel[j]), chromend = end(reducedDel[j]),
      y = subplotHeight, height = currentHeight-subplotHeight
    )
  } else {
    annoHighlight(
      plot=nhr46Plot,
      params=nhr46params,
      chrom = as.character(seqnames(reducedDel[j])),
      chromstart = start(reducedDel[j]), chromend = end(reducedDel[j]),
      y = subplotHeight, height = currentHeight-subplotHeight
    )
  }
}


annoGenomeLabel(
  plot = skn1Plot,
  params=skn1params,
  y = currentHeight, height = textHeight,
)

annoGenomeLabel(
  plot = nhr46params,
  params=nhr46params,
  y = currentHeight, height = textHeight,
)

plotRect(x=0.5,y=subplotHeight-hicHeight+2,
         width=10,height=currentHeight-subplotHeight + hicHeight -2,
         params=skn1params)
#y=11.3 height=13
plotRect(x=11,y=subplotHeight-hicHeight+2,
         width=6.5,height=currentHeight-subplotHeight + hicHeight -2,
         params=nhr46params)

dev.off()

