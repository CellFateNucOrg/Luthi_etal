library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(plyranges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggtext)
library(cowplot)
library(ggbeeswarm)

theme_set(
  theme_bw(base_size=9)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=9),
          plot.title=element_text(size=9),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          strip.text=element_text(size=9),
          legend.title=element_text(size=9)
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


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125b.rds"))
seqlevels(fountains)<-seqlevels(Celegans)
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")


# control bins
nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=2000,fix="center")
rtracklayer::export(nonFount,paste0(fountainsDir,"/controlsForFountains_base0_uncorrected_20240125.bed"))

# enhancers
daughertyL3<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))


metadata<-readRDS("/Users/semple/Documents/MeisterLab/papers/Isiaka_etal/wbGeneGR_WS285.rds")
seqlevelsStyle(metadata)<-"UCSC"

operons<-import("/Users/semple/Documents/MeisterLab/papers/Moushumi1/tracks/operon.bed")
seqlevelsStyle(operons)<-"UCSC"
seqlevels(operons)<-seqlevels(Celegans)
sum(countOverlaps(fountains,operons))
sum(countOverlaps(nonFount,operons))

df<-data.frame(class=table(metadata$class))
colnames(df)<-c("class","count")
df$fount<-0
df$cont<-0
for(c in names(table(metadata$class))){
  print(c)
  df$fount[df$class==c]<-sum(countOverlaps(fountains,metadata[metadata$class==c]))
  df[df$class==c,"cont"]<-sum(countOverlaps(nonFount,metadata[metadata$class==c]))
}
df

ol<-findOverlaps(fountains,metadata[metadata$class=="ncRNA_gene"])
hist(width(metadata[metadata$class=="ncRNA_gene"][subjectHits(ol)]),breaks=100)
metadata[metadata$class=="piRNA_gene"]


operons<-import("/Users/semple/Documents/MeisterLab/papers/Moushumi1/tracks/operon.bed")
seqlevelsStyle(operons)<-"UCSC"
seqlevels(operons)<-seqlevels(Celegans)
sum(countOverlaps(fountains,operons)) #163
sum(countOverlaps(nonFount,operons)) #185

