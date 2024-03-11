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
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=0.4,size=3,colour="blue") +
  xlab("Distance from fountain bins (kb)") + ylab("log2(baseMean +1)") +
  ggtitle("Log2 fold change coh1cs/TEV control by distance from fountain") +
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
  xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+1)")+
  ggtitle("Log2 TPM in TEV control animals by distance from fountain") +
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
  xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+1)")+
  ggtitle("Log2 TPM in COH-1 cleaved animals by distance from fountain") +
  coord_cartesian(ylim=c(0,13))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
                     hide.ns = T, color="purple",x="group2",y=11)+
  geom_hline(yintercept=mean(df[selected,"log2TPMcoh1"]),linetype="dashed")
p3
stat_df


p<-ggarrange(p1,p2,p3,nrow=3)
p<-annotate_figure(p, top = text_grob("Isiaka et al., Supl. Figure", size = 14))
ggsave(filename=paste0(finalFigDir,"/supplFig_fountainExpression_2kb_bins.pdf"),
       p, device="pdf",height=29,width=19, units="cm")

