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
fountainFigDir=paste0(projectDir,"/fountainFigures")

fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
seqlevels(fountains)<-seqlevels(Celegans)
seqinfo(fountains)<-seqinfo(Celegans)
i=500
resolution=2000
region=3e5
pageCreate(width = 18, height = 27, params=params, showGuides=FALSE)

drawGenomicRegion<-function(fountains,fountainIndex,resolution=resolution,
                            regionSize=region,currentHeight=0){
  gr<-trim(resize(fountains[i],width=max(regions)+resolution,fix="center"))
  print(paste0("Fountain: ", as.character(fountains[i])))
  print(paste0("Gathering data for: ",as.character(gr)))
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

  ## plotgardener framework
  mat<-data.frame(chr=start(anchors(div_contacts1)$first),chr=start(anchors(div_contacts1)$second),count=scores(div_contacts1)$balanced.l2fc)

  #change gr to fit hic matrix
  grFix<-GRanges(seqnames=seqnames(gr),ranges=IRanges(start=min(mat$chr),end=max(mat$chr.1)))
  seqlevels(grFix)<-seqlevels(Celegans)
  seqinfo(grFix)<-seqinfo(Celegans)

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


  chromDomains<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromDomains_L3_Evans2016_ce11_rgb.bed")





  ############ plots -----
  grFix<-trim(resize(fountains[i],width=region+resolution,fix="center"))
  print(paste0("Plotting: ",as.character(grFix)))
  if(!dir.exists(paste0(fountainFigDir,"/reg",region/1e3,"kb_res",resolution/1e3,"kb/"))){
    dir.create(paste0(fountainFigDir,"/reg",region/1e3,"kb_res",resolution/1e3,"kb/"),
               recursive=T)
  }

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
  hicHeight=5
  transcriptHeight=3
  chipBin=200

  hicPlot<-plotHicTriangle(
    data = mat,
    params=params,
    palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
    zrange=c(-1,1),
    x = 0.5, y = currentHeight, width =17, height=hicHeight,
    colorTrans="linear"
  )

  annoHeatmapLegend(
    plot = hicPlot,
    params=params,
    x = 17, y = currentHeight+0.2, width = 0.3, height = 2,
    fontcolor="black"
  )

  plotText(label = "HiC\nlog2(COH-1/TEVonly)", fontcolor = "black",
           y = currentHeight+0.2, params=params)

  currentHeight<-currentHeight+ hicHeight  + 0.2

  ## fountains -----
  p<-plotRanges(fountains,
                params=params,
                y = currentHeight, height=rangesHeight,
                boxHeight=boxHeight,
                fill="red")
  plotText(label = "Fountain tip", fontcolor = "red",
           y = currentHeight+0.1, params=params)
  currentHeight<-currentHeight+rangesHeight+0.1


  ## enahncers -----
  p<-plotRanges(daughertyL3,
                params=params,
                y = currentHeight, height=rangesHeight,
                boxHeight=boxHeight,
                spaceHeight=0,spaceWidth=0,
                fill="violet")
  plotText(label = "Enhancers(Daugherty)", fontcolor = "black",
           y = currentHeight+0.1, params=params)
  currentHeight<-currentHeight+rangesHeight + 0.1

  p<-plotRanges(JaenesL3,
                params=params,
                y = currentHeight, height=rangesHeight,
                spaceHeight=0,spaceWidth=0,
                boxHeight=boxHeight,
                fill="violet")
  plotText(label = "Enhancers(Jaenes)", fontcolor = "black",
           y = currentHeight+0.1, params=params)
  currentHeight<-currentHeight+rangesHeight + 0.1

  ## ATAC and chip seq----
  p<-plotSignal(coh1,
                params=params,
                binSize=chipBin,
                fill="red2",
                y = currentHeight, height=signalHeight,
                linecolor="red2")
  plotText(label = "COH-1 ChIP", fontcolor = "red2",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1

  p<-plotSignal(daughertyATAC,
                params=params,
                binSize=chipBin,
                fill="violet",
                y = currentHeight, height=signalHeight,
                linecolor="violet")
  plotText(label = "ATAC(Daugherty)", fontcolor = "violet",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1


  p<-plotSignal(jaenesATAC,
                params=params,
                binSize=chipBin,
                fill="violet",
                y = currentHeight, height=signalHeight,
                linecolor="violet")
  plotText(label = "ATAC(Jaenes)", fontcolor = "violet",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1


  p<-plotSignal(H3K27ac,
                params=params,
                binSize=chipBin,
                fill="darkgreen",
                y = currentHeight, height=signalHeight,
                linecolor="darkgreen",
                range=c(0,max(H3K27ac$score)),negData=T)
  plotText(label = "H3K27Ac ChIP", fontcolor = "darkgreen",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1

  p<-plotSignal(H3K27me3,
                params=params,
                binSize=chipBin,
                fill="darkgrey",
                y = currentHeight, height=signalHeight,
                linecolor="darkgrey",
                range=c(0,max(H3K27me3$score)))
  plotText(label = "H3K27me3 ChIP", fontcolor = "darkgrey",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1


  p<-plotSignal(top1,
                params=params,
                binSize=chipBin,
                fill="orange2",
                y = currentHeight, height=signalHeight,
                linecolor="orange2",
                negData=T)
  plotText(label = "TOP-1 ChIP", fontcolor = "orange2",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1

  p<-plotSignal(top2,
                params=params,
                binSize=chipBin,
                fill="orange4",
                y = currentHeight, height=signalHeight,
                linecolor="orange4",
                negData=T)
  plotText(label = "TOP-2 ChIP", fontcolor = "orange4",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + 0.1

  ## chromHMM -----
  plotText(label = "ChromatinDomains (Evans) ", fontcolor = "black",
           y = currentHeight, params=params)
  p<-plotRanges(chromDomains,
                params=params,
                y = currentHeight+textHeight, height=rangesHeight,
                boxHeight=boxHeight, fill=chromDomains$itemRgb)
  currentHeight<-currentHeight + textHeight+rangesHeight + 0.1

  p<-plotGenes(params=params,
               y=currentHeight, height=signalHeight,
               fontsize=1)
  currentHeight<-currentHeight + signalHeight + 0.1

  fountLoc<-subsetByOverlaps(fountains,grFix)

  for(j in 1:length(fountLoc)){
    annoHighlight(
      plot=p,
      params=params,
      chrom = as.character(seqnames(fountLoc[j])),
      chromstart = start(fountLoc[j]), chromend = end(fountLoc[j]),
      y = hicHeight, height = currentHeight-hicHeight
    )
  }

  annoGenomeLabel(
    plot = p,
    params=params,
    y = currentHeight, height = textHeight,
  )
}
