library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(plyranges)
library(rtracklayer)
library(dplyr)
library(ggtext)
library(cowplot)

library(plotgardener)
library(HiCExperiment)
library(HiContacts)
library(RColorBrewer)
library(TxDb.Celegans.UCSC.ce11.refGene)



theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9)
    )
)


projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseq_DGE")
rnaTxSeqDir=paste0(projectDir,"/RNAseq_DTE")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

source(paste0(projectDir,"/functions.R"))
source(paste0(projectDir,"/functions_fountainPlots.R"))


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
# fountains$AsymmetryTercile<-paste0("Asym.Q",fountains$AsymmetryTercile)
# fountains$AsymmetryQuintile<-paste0("Asym.Q",fountains$AsymmetryQuintile)





######################################
## chromatin domains ---------
######################################

domains<-import.bed(paste0(publicDataDir,"/chromDomains_L3_Evans2016_ce11.bed"))
seqlevels(domains)<-seqlevels(Celegans)
domains$name<-factor(domains$name)
domains$score<-factor(as.numeric(domains$name))

domainClrs<-c("#cc0000","#bcbcbc","#5b5b5b")

# getDomainOLtable<-function(gr,domains){
#   ol<-data.frame(findOverlaps(gr,domains,ignore.strand=T,minoverlap=10))
#   gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
#   ol$name<-gr$name[ol$queryHits]
#   ol$domain<-factor(domains$score[ol$subjectHits],levels=1:3)
#   ol$width<-width(domains)[ol$subjectHits]
#   df<-ol%>%dplyr::group_by(domain) %>%
#     dplyr::summarise(domainFrequency=dplyr::n(),domainWidth=sum(width))
#   alldomains<-data.frame(domain=factor(1:3, levels=1:3))
#   df<-left_join(alldomains,df)
#   df[is.na(df)]<-0
#   return(df)
# }


getDomainOLtable<-function(gr,domains){
  ol<-data.frame(findOverlaps(gr,domains,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$domain<-factor(domains$name[ol$subjectHits],levels=c("active","border","regulated"))
  ol$width<-width(domains)[ol$subjectHits]
  df<-ol%>%dplyr::group_by(domain) %>%
    dplyr::summarise(domainFrequency=dplyr::n(),domainWidth=sum(width))
  df[is.na(df)]<-0
  return(df)
}

tiles<-tileGenome(seqlengths(Celegans)[1:6],tilewidth=2000,cut.last.tile.in.chrom = T)
tiles<-getUpDownstreamFountainData(tiles,fountains)

tiles<-makeGRangesFromDataFrame(tiles,keep.extra.columns = T)

tiles$binnedDistance<-binByDistance1(tiles$distanceToFountain,maxDist=20000,binSize =2000,fullTile=T)
tiles<-flipDownstreamOrientation(tiles)
table(tiles$binnedDistance1,tiles$fountainLocation)
selected<-levels(tiles$binnedDistance1)[-c(grep(">",levels(tiles$binnedDistance1)),grep("-0",levels(tiles$binnedDistance1)))]

listdf<-list()
for(l in levels(tiles$fountainLocation)){
    for(d in selected){
      df1<-getDomainOLtable(gr=tiles[tiles$binnedDistance1==d],domains)
      df1$fountainLocation<-l
      df1$binnedDistance1<-d
      listdf[[paste0(l,"_",d)]]<-df1
    }
}

df<-do.call(rbind,listdf)


df$binnedDistance1<-factor(df$binnedDistance1,levels=selected)

bin0<-ceiling(length(levels(df$binnedDistance1))/2)

boxdf<-df %>%
  dplyr::group_by(binnedDistance1) %>%
  mutate(domainWidth=sum(domainWidth),domainFrequency=sum(domainFrequency)) %>%
  filter(binnedDistance1=="0") %>% dplyr::select(-domain,-fountainLocation) %>% distinct()

p1<-ggplot(df,aes(x=binnedDistance1,y=domainFrequency/1e3,fill=domain)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_manual(values=domainClrs) +
  xlab("Binned distance (kb)")+
  theme(legend.position = "bottom",legend.title=element_blank(),
        legend.key.size = unit(0.3, 'cm'),legend.box.spacing = unit(0, "pt")) +
  #ggtitle(paste0("Chromatin domain frequency around autosomal fountain tips")) +
  ylab("Bin domain overlap frequency (x10<sup>3</sup>)")
p1<-p1 + geom_bar(data=boxdf,aes(x=binnedDistance1, y=domainFrequency/1e3,fill=NULL),stat="identity",
                  alpha=0,color="black",linewidth=0.5)

p1


# genome average
domains$width=width(domains)
genomeWide<-data.frame(domains) %>% group_by(name, score) %>% summarise(domainWidth=sum(width),
                                                                        domainFreq=n())

p2<-ggplot(genomeWide,aes(x=1,y=domainWidth/1e7,fill=name)) +
  geom_bar(position="stack",stat="identity")  +
  scale_fill_manual(values=domainClrs) +
  ylab("Genome-wide domain width (x10<sup>7</sup>)") + ggtitle("")+
  theme(legend.position="none",title=element_blank(),axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,"cm"),axis.text.x=element_blank())
p2



ol<-findOverlaps(fountains,domains[domains$name=="active"])

fountains$active<-countOverlaps(fountains,domains[domains$name=="active"])
fountains$border<-countOverlaps(fountains,domains[domains$name=="border"])
fountains$regulated<-countOverlaps(fountains,domains[domains$name=="regulated"])
tb<-tidyr::pivot_longer(data.frame(fountains), cols=c("active","border","regulated"),
                    names_to="type",values_to="counts")

stat_box_data <- function(y,upper_limit=0) {
  return(
    data.frame(
      y = upper_limit,
      label = paste(length(y))
    )
  )
}

# Statistical test
stat.test <- tb[tb$counts>0,] %>%
  wilcox_test(fountain.score ~ type) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "type")
stat.test$p.format <- p_format(
  stat.test$p.adj, accuracy = 1e-7,
  leading.zero = T
)


p3<-ggplot(tb[tb$counts>0,],aes(x=type,y=fountain.score*1e5)) +
  geom_boxplot(aes(fill=type),notch=T,outlier.shape=NA) +
  coord_cartesian(ylim=c(0,3)) +
  scale_fill_manual(values=domainClrs) +
  theme(legend.position = "none", axis.title.x=element_blank(),
        legend.title=element_blank(),title=element_blank())+
  ylab("fountain score (x10<sup>-5</sup>)")+
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.9,
    size = 3
  )

p3<-p3 + stat_pvalue_manual(data=stat.test, label = "p.format",y.position=c(3,2.8,2.6),
                       tip.length=0.005,color="purple",size=3)

p<-cowplot::plot_grid(p1,p2,p3,nrow=1,ncol=3,rel_widths=c(0.9,0.12,0.3),labels=c("a","","b"),
                      align="h")

p<-annotate_figure(p, top = text_grob("Lüthi et al., Figure S4", size = 14))
#ggsave(paste0(finalFigDir,"/supplFig_EvansDomains_2kb_bins_freq.pdf"), p, device="pdf",
#       width=19,height=10, units="cm")

p


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
seqlevels(fountains)<-seqlevels(Celegans)
seqinfo(fountains)<-seqinfo(Celegans)
resolution=2000
region=3e5


drawGenomicRegion<-function(fountains,fountainIndex,resolution=resolution,
                            regionSize=region,initialHeight=0){
  currentHeight<-initialHeight
  i=fountainIndex
  gr<-trim(resize(fountains[i],width=region+resolution,fix="center"))
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
  #GSE50324

  H3K27ac<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/GSM624432_WA30634849_H3K27AC_N2_L3_1_ce11.bw",
                  selection=grFix)
  #GSM624432

  #H3K27me3<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/GSM562734_HK00013_H3K27ME31E7_N2_L3_1_ce11.bw",
  #                 selection=grFix)

  H3K27me3<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_realigned/modEncode_histone/bigwig_fc/H3K27me3_UP07449_rep.pooled_x_ctl.pooled.fc.signal.bigwig",
                   selection=grFix)
#GSM1206310

  top1<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP1-GFP_AM05_L2-L3_ce11_avr.bw",
               selection=grFix)
 # average of GSM5686806 and GSM5686807
 #
  top2<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Morao-Ercan_MolCell2022_GSE188851/averaged/TOP2-GFP_MDX53_L2-L3_ce11_avr.bw",
               selection=grFix)
  #average of GSM5686812 and GSM5686813

  chromDomains<-import("/Volumes/external.data/MeisterLab/publicData/ChIPseq_liftedOver/Evans2016_PNAS_L3/chromDomains_L3_Evans2016_ce11_rgb.bed")





  ############ plots -----
  grFix<-trim(resize(fountains[i],width=region+resolution,fix="center"))
  print(paste0("Plotting: ",as.character(grFix)))
  if(!dir.exists(paste0(fountainsDir,"/reg",region/1e3,"kb_res",resolution/1e3,"kb/"))){
    dir.create(paste0(fountainsDir,"/reg",region/1e3,"kb_res",resolution/1e3,"kb/"),
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
  spacerHeight=0.1

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

  plotText(label = "log2(COH-1/TEVonly)\nHiC", fontcolor = "black",
           y = currentHeight+0.5, params=params)

  currentHeight<-currentHeight+ hicHeight + spacerHeight

  ## fountains -----
  p<-plotRanges(fountains,
                params=params,
                y = currentHeight, height=rangesHeight,
                boxHeight=boxHeight,
                fill="red")
  plotText(label = "Fountain tip", fontcolor = "red",
           y = currentHeight+0.1, params=params)
  currentHeight<-currentHeight+rangesHeight + spacerHeight


  ## enahncers -----
  p<-plotRanges(daughertyL3,
                params=params,
                y = currentHeight, height=rangesHeight,
                boxHeight=boxHeight,
                spaceHeight=0,spaceWidth=0,
                fill="violet")
  plotText(label = "Enhancers(Daugherty)", fontcolor = "black",
           y = currentHeight+0.1, params=params)
  currentHeight<-currentHeight+rangesHeight + spacerHeight

  p<-plotRanges(JaenesL3,
                params=params,
                y = currentHeight, height=rangesHeight,
                spaceHeight=0,spaceWidth=0,
                boxHeight=boxHeight,
                fill="violet")
  plotText(label = "Enhancers(Jaenes)", fontcolor = "black",
           y = currentHeight+0.1, params=params)
  currentHeight<-currentHeight+rangesHeight + spacerHeight

  ## ATAC and chip seq----
  p<-plotSignal(coh1,
                params=params,
                binSize=chipBin,
                fill="red2",
                y = currentHeight, height=signalHeight,
                linecolor="red2")
  plotText(label = "COH-1 ChIP", fontcolor = "red2",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight

  p<-plotSignal(daughertyATAC,
                params=params,
                binSize=chipBin,
                fill="violet",
                y = currentHeight, height=signalHeight,
                linecolor="violet")
  plotText(label = "ATAC(Daugherty)", fontcolor = "violet",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight


  p<-plotSignal(jaenesATAC,
                params=params,
                binSize=chipBin,
                fill="violet",
                y = currentHeight, height=signalHeight,
                linecolor="violet")
  plotText(label = "ATAC(Jaenes)", fontcolor = "violet",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight


  p<-plotSignal(H3K27ac,
                params=params,
                binSize=chipBin,
                fill="darkgreen",
                y = currentHeight, height=signalHeight,
                linecolor="darkgreen",
                range=c(0,max(H3K27ac$score)),negData=T)
  plotText(label = "H3K27Ac ChIP", fontcolor = "darkgreen",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight

  p<-plotSignal(H3K27me3,
                params=params,
                binSize=chipBin,
                fill="darkgrey",
                y = currentHeight, height=signalHeight,
                linecolor="darkgrey",
                range=c(0,max(H3K27me3$score)))
  plotText(label = "H3K27me3 ChIP", fontcolor = "darkgrey",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight


  p<-plotSignal(top1,
                params=params,
                binSize=chipBin,
                fill="orange2",
                y = currentHeight, height=signalHeight,
                linecolor="orange2",
                negData=T)
  plotText(label = "TOP-1 ChIP", fontcolor = "orange2",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight

  p<-plotSignal(top2,
                params=params,
                binSize=chipBin,
                fill="orange4",
                y = currentHeight, height=signalHeight,
                linecolor="orange4",
                negData=T)
  plotText(label = "TOP-2 ChIP", fontcolor = "orange4",
           y = currentHeight, params=params)
  currentHeight<-currentHeight + signalHeight + spacerHeight

  ## chromHMM -----
  plotText(label = "ChromatinDomains (Evans) ", fontcolor = "black",
           y = currentHeight, params=params)
  p<-plotRanges(chromDomains,
                params=params,
                y = currentHeight+textHeight, height=rangesHeight,
                boxHeight=boxHeight, fill=chromDomains$itemRgb)
  currentHeight<-currentHeight + textHeight+rangesHeight + spacerHeight

  p<-plotGenes(params=params,
               y=currentHeight, height=signalHeight,
               fontsize=0)
  currentHeight<-currentHeight + signalHeight + spacerHeight

  fountLoc<-subsetByOverlaps(fountains,grFix)

  for(j in 1:length(fountLoc)){
    annoHighlight(
      plot=p,
      params=params,
      chrom = as.character(seqnames(fountLoc[j])),
      chromstart = start(fountLoc[j]), chromend = end(fountLoc[j]),
      y = initialHeight+hicHeight, height = currentHeight-hicHeight-initialHeight
    )
  }

  annoGenomeLabel(
    plot = p,
    params=params, scale="bp",
    y = currentHeight, height = textHeight,
  )
}

pdf(paste0(finalFigDir,"/FigS4_EvansDomains_2kb_bins_freq.pdf"),
       width=11,height=19,paper="a4")

pageCreate(width = 18, height = 25, params=params, showGuides=FALSE)
plotGG(plot=p,x=0.5,y=0,width=17,height=10,default.units="cm")
plotText(label = "c", fontsize = 16, fontface = "bold",
  x = 0.5, y = 10, just = "center", default.units = "cm")
drawGenomicRegion(fountains,fountainIndex=249,resolution=resolution,
                            regionSize=region,initialHeight=10)
dev.off()

table(countOverlaps(domains[domains$name=="active"],fountains))
