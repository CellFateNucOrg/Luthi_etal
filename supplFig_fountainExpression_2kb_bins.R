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


theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
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

rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")




## log2 FC ------
salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

salmongr<-sort(salmongr)
salmongr$geneLength<-width(salmongr)
salmongr<-resize(salmongr,width=1,fix="start")
df<-getUpDownstreamFountainData(salmongr,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=20000,binSize=2000)
table(df$binnedDistance,df$AsymmetryQuintile)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)

ignoreBins<-c(">-20","-0",">20")

selected<-!df$binnedDistance1 %in% ignoreBins
countLabels_a<-df[selected,]%>% dplyr::group_by(binnedDistance1) %>% dplyr::summarise(count=dplyr::n())

# calculate stats
stat_df<-df %>%
  rstatix::t_test(log2FoldChange~binnedDistance1,ref.group="all")  %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group2 %in% ignoreBins) %>%
  add_x_position(x="binnedDistance1")

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p1<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2FoldChange),
               fill=c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
               outlier.color=NA,notch=T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1), strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=0.4,size=3,colour="blue") +
  xlab("Distance from fountain bins (kb)") + ylab("log<sub>2</sub>FC") +
  ggtitle("Log<sub>2</sub> fold change coh1cs/TEV control by distance from fountain") +
  stat_pvalue_manual(stat_df,label = "p.adj.signif", remove.bracket=T,hide.ns = T,
                     color="purple",x="group2",y=0.35,size=7) + coord_cartesian(ylim=c(-0.25,0.4))+
  geom_hline(yintercept=0,linetype="dashed")
p1


## 366 TPM - 10 kb zoom----
tpm366<-import(paste0(rnaSeqDir,"/tracks/PMW366_TPM_avr.bed"))
tpm366$length<-width(tpm366)
tpm366tss<-resize(tpm366,width=1,fix="start")
df<-getUpDownstreamFountainData(tpm366tss,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=20000,binSize=2000)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)

df$log2TPMtev<-log2(df$score+1)
table(df$binnedDistance,df$AsymmetryTercile)
table(df$binnedDistance1)


ignoreBins<-c(">-20","-0",">20")

selected=!df$binnedDistance1 %in% ignoreBins & df$log2TPMtev >0
countLabels_a<-df[selected,]%>%
  dplyr::group_by(binnedDistance1,.drop=T) %>%
  dplyr::summarise(count=dplyr::n())

# calculate stats
stat_df<-df %>%
  rstatix::t_test(log2TPMtev~binnedDistance1,ref.group="all")  %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group2 %in% ignoreBins) %>%
  add_xy_position(x="binnedDistance1")

stat_df

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p2<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2TPMtev),
               fill=c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
               outlier.color=NA) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=12,size=3,colour="blue",
            angle=0) +
  xlab("Distance from fountain bins (kb)") + ylab("Log<sub>2</sub>(TPM+1)")+
  ggtitle("Log<sub>2</sub> TPM in TEV control animals by distance from fountain") +
  coord_cartesian(ylim=c(0,13))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=T,
                     hide.ns = T, color="purple",x="group2",y=11) +
  geom_hline(yintercept=mean(df[selected,"log2TPMtev"]),linetype="dashed")
p2
stat_df$p.adj





## coh-1 TPM - 10 kb zoom----
tpmCOH1<-import(paste0(rnaSeqDir,"/tracks/PMW828_TPM_avr.bed"))

tpmCOH1$length<-width(tpmCOH1)
tpmCOH1tss<-resize(tpmCOH1,width=1,fix="start")
df<-getUpDownstreamFountainData(tpmCOH1tss,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=20000,binSize=2000)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)


table(df$binnedDistance,df$AsymmetryTercile)

df$log2TPMcoh1<-log2(df$score+1)
ignoreBins<-c(">-20","-0",">20")


selected=!df$binnedDistance1 %in% ignoreBins & df$log2TPMcoh1 >0
countLabels_a<-df[selected,]%>%
  dplyr::group_by(binnedDistance1,.drop=T) %>%
  dplyr::summarise(count=dplyr::n())

# calculate stats
stat_df<-df[df$log2TPMcoh1 >0,] %>%
  rstatix::t_test(log2TPMcoh1~binnedDistance1,ref.group="all")  %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group2 %in% ignoreBins) %>%
  add_xy_position(x="binnedDistance1")

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p3<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2TPMcoh1),
               fill=c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
               outlier.color=NA) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=12,size=3,colour="blue",
            angle=0) +
  xlab("Distance from fountain bins (kb)") + ylab("Log<sub>2</sub>(TPM+1)")+
  ggtitle("Log<sub>2</sub> TPM in COH-1 cleaved animals by distance from fountain") +
  coord_cartesian(ylim=c(0,13))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
                     hide.ns = T, color="purple",x="group2",y=11)+
  geom_hline(yintercept=mean(df[selected,"log2TPMcoh1"]),linetype="dashed")
p3
stat_df



## number of enhancers per bin vs LFC
# enhancers
daughertyL3<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))

binSize=6000
fountains$fountVcont<-"fountain tip"
#nonFount$fountVcont<-"control"

allRegions<-fountains
allRegions<-resize(allRegions,width=binSize,fix="center")

allRegions$ActiveCount<-countOverlaps(allRegions,daughertyL3[daughertyL3$L3_chromHMMState=="L3_activeEnhancer"],ignore.strand=T,minoverlap=1)

salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
salmonGR<-GRanges(salmon)
ol<-findOverlaps(resize(salmonGR,width=1,fix="start"),allRegions,ignore.strand=T)

salmonGR$ActiveCount<-0
salmonGR$ActiveCount[queryHits(ol)]<-allRegions$ActiveCount[subjectHits(ol)]

df<-data.frame(salmonGR)
df$ActiveCount<-factor(df$ActiveCount)

#ggplot(data=df[df$ActiveCount!=0,],aes(x=ActiveCount,y=log2FoldChange,fill=ActiveCount)) + geom_beeswarm() +
#  geom_hline(yintercept=0) + scale_fill_grey(start=1,end=0.4)


#' Function for adding count data to plots
count_data <- function (y,ymax=5){
  df <- data.frame(y = ymax, label = length(y))
  return(df)
}


p4<-ggplot(data=df,aes(x=ActiveCount,y=log2FoldChange)) +
  geom_jitter(width=0.2,size=1,aes(color=ActiveCount)) +
  geom_boxplot(outlier.shape=NA,fill=NA,width=0.3) +
  geom_hline(yintercept=0,color="red",linetype="dashed") +
  scale_color_grey(start=1,end=0) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  xlab(paste0("Number of active enhancers per ",binSize/1000,"kb fountain tip bin")) +
  ylab(label="Log<sub>2</sub>FC") +
  stat_summary(fun.data = count_data,fun.args=c(ymax=0.5), geom = "text",
               position = position_dodge(1), size=3, colour="blue")+
  theme(legend.position="none") +
  ggtitle("L3 enhancers (Daugherty et al. 2017)")
p4




## number of enhancers per bin vs LFC

# enhancers
JaenesL<-readRDS(paste0(publicDataDir,"/Jaenes2018_enhancers_ce11_stages_chromHMM.rds"))
JaenesL<-JaenesL%>% filter(topState_L3_chromHMM %in% c("Active enhancer","H3K27me3 repressed","Repressed enhancer"))
length(JaenesL)

binSize=6000
fountains$fountVcont<-"fountain tip"
#nonFount$fountVcont<-"control"

allRegions<-fountains
allRegions<-resize(allRegions,width=binSize,fix="center")

allRegions$ActiveCount<-countOverlaps(allRegions,JaenesL[JaenesL$topState_L3_chromHMM=="Active enhancer"],ignore.strand=T,minoverlap=1)

salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
salmonGR<-GRanges(salmon)
ol<-findOverlaps(resize(salmonGR,width=1,fix="start"),allRegions,ignore.strand=T)

salmonGR$ActiveCount<-0
salmonGR$ActiveCount[queryHits(ol)]<-allRegions$ActiveCount[subjectHits(ol)]

df<-data.frame(salmonGR)
df$ActiveCount<-factor(df$ActiveCount)


#' Function for adding count data to plots
count_data <- function (y,ymax=5){
  df <- data.frame(y = ymax, label = length(y))
  return(df)
}


p5<-ggplot(data=df,aes(x=ActiveCount,y=log2FoldChange)) +
  geom_jitter(width=0.2,size=1,aes(color=ActiveCount)) +
  geom_boxplot(outlier.shape=NA,fill=NA,width=0.3) +
  geom_hline(yintercept=0,color="red",linetype="dashed") +
  scale_color_grey(start=1,end=0) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  xlab(paste0("Number of active enhancers per ",binSize/1000,"kb fountain tip bin")) +
  ylab(label="Log<sub>2</sub>FC") +
  stat_summary(fun.data = count_data,fun.args=c(ymax=0.5), geom = "text",
               position = position_dodge(1), size=3, colour="blue")+
  theme(legend.position="none") +
  ggtitle("L3 enhancers (Jaenes et al. 2018)")
p5








p<-ggarrange(p1,p2,p3,ggarrange(p4,p5,nrow=1,ncol=2,labels=c("d ","e ")),nrow=4,
             labels=c("a ","b ","c "))
p<-annotate_figure(p, top = text_grob("Lüthi et al., Figure S7", size = 12))
ggsave(filename=paste0(finalFigDir,"/supplFig_fountainExpression_2kb_bins.pdf"),
       p, device="pdf",height=29,width=19, units="cm")

