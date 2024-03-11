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

theme_set(
  theme_bw(base_size=9)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=9),
          plot.title=element_text(size=9),
          axis.title=element_text(size=9),
          strip.text=element_text(size=9)
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

# control bins
nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=2000,fix="center")


#########################-
# Daugherty ATAC enhancers -----
#########################-


# enhancers
daughertyL3<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))


####### ecdf enhancer distance to fountains
daughertyL3$distanceToFount<-mcols(distanceToNearest(daughertyL3,fountains,ignore.strand=T))$distance

options(tibble.width=Inf)
dd1<-data.frame(daughertyL3) %>%
  dplyr::group_by(L3_chromHMMState) %>%
  dplyr::mutate(ecd=ecdf(distanceToFount)(distanceToFount))

dd1$type<-factor(dd1$L3_chromHMMState, levels=c("L3_activeEnhancer" ,"L3_repressedEnhancer","L3_H3K27me3Repressed"))
levels(dd1$type)<-c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer")

dd1$type_count<-dd1$type
levels(dd1$type_count)<-paste0(c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer"),
                               paste0(" (n=",dd2$count,")"))

p1<-ggplot(dd1, aes(x=distanceToFount/1000,y=ecd,color=type_count)) +
  geom_line(linewidth=0.9)+
  scale_color_manual(values=c("red","blue","darkgreen"))+
  xlab("Distance to fountain tip (kb)")+ylab("Cumulative distribution of distance")+
  coord_cartesian(xlim=c(1,150000/1000)) +
  geom_hline(yintercept=1,colour="darkgrey") +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.5)) +
  ggtitle("Daugherty et al. (2017) L3 enhancers")

# to get fraction between 0.15 and 0.5 lfc
dd2<-dd1 %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(ecd10kb=ecdf(distanceToFount)(c(10000)),
                   ecd50kb=ecdf(distanceToFount)(c(50000)),
                   bracket1=distanceToFount[which.min(abs(ecd-0.6))],
                   bracket2=distanceToFount[which.min(abs(ecd-0.63))],
                   count=n())
# add 10kb quantile
p1<-p1+ geom_segment(x=10,y=dd2$ecd10kb[dd2$type=="Active enhancer"],xend=10,
                     yend=-Inf,linetype="dashed",colour="grey") +
  geom_segment(x=-Inf,y=dd2$ecd10kb[dd2$type=="Active enhancer"],
               xend=10,yend=dd2$ecd10kb[dd2$type=="Active enhancer"],
               linetype="dashed",colour="grey") +
  annotate(geom="text",x=11,y=0,label="10kb",color="darkgrey",hjust=0,size=3)+
  annotate(geom="text",x=0,y=dd2$ecd10kb[dd2$type=="Active enhancer"]*1.02,
           label=paste0(round((dd2$ecd10kb[dd2$type=="Active enhancer"])*100,0),"%"),
           color="darkgrey",vjust=0,size=3)
# add pvalues
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="Repressed enhancer"],alternative="less")
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")
ks.test(dd1$ecd[dd1$type=="Repressed enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")

p1<-p1+ geom_bracket(xmin=as.numeric(dd2[1,"bracket1"]/1000),
                     xmax=as.numeric(dd2[2,"bracket1"]/1000),
                     y.position=0.6, label="****",inherit.aes=F, bracket.nudge.y=0.02,
                     tip.length=0.02,vjust=0.75,label.size=3) +
  geom_bracket(xmin=as.numeric(dd2[2,"bracket2"]/1000),
               xmax=as.numeric(dd2[3,"bracket2"]/1000),
               y.position=0.63, label="ns",inherit.aes=F, bracket.nudge.y=0.02,
               tip.length=0.02,vjust=0.2,label.size=3)



# Number of enhancers per bin -------
binSize=6000
fountains$fountVcont<-"fountain tip"
nonFount$fountVcont<-"control"

allRegions<-c(fountains,nonFount)
allRegions<-resize(allRegions,width=binSize,fix="center")

allRegions$ActiveCount<-countOverlaps(allRegions,daughertyL3[daughertyL3$L3_chromHMMState=="L3_activeEnhancer"],ignore.strand=T,minoverlap=1)
allRegions$RepressedCount<-countOverlaps(allRegions,daughertyL3[daughertyL3$L3_chromHMMState=="L3_repressedEnhancer"],ignore.strand=T,minoverlap=1)
allRegions$H3K27meRepressedCount<-countOverlaps(allRegions,daughertyL3[daughertyL3$L3_chromHMMState=="L3_H3K27me3Repressed"],ignore.strand=T,minoverlap=1)


df<-data.frame(allRegions)

df1<-pivot_longer(df,cols=c("ActiveCount", "RepressedCount", "H3K27meRepressedCount"),
             names_to="type",values_to="count")

df1$type<-gsub("Count$","",df1$type)
df1$type<-factor(df1$type,levels=c("Active", "Repressed","H3K27meRepressed"))
levels(df1$type)<-c("Active", "Repressed","H3K27me3")

p2<-ggplot(df1,aes(x=count,fill=type)) +
  geom_bar(position=position_dodge(preserve = "single")) +
  facet_wrap(.~fountVcont) +
  scale_fill_manual(values=c("red","blue","darkgreen")) +
  xlab(paste0("Number of enhancers per ",binSize/1000,"kb bin"))+
  ylab(paste0("Number of ",binSize/1000,"kb bins"))+
  theme(legend.title=element_blank(), legend.position=c(0.8,0.7),
        legend.key.size = unit(0.3, 'cm'))+
  ggtitle("Daugherty et al. (2017) L3 enhancers")




#########################-
# Jaenes ATAC enhancers -----
#########################-


# enhancers
JaenesL<-readRDS(paste0(publicDataDir,"/Jaenes2018_enhancers_ce11_stages_chromHMM.rds"))
JaenesL<-JaenesL%>% filter(maxStage!="emb",maxStage!="l4",maxStage!="ya", topState_L3_chromHMM %in% c("Active enhancer","H3K27me3 repressed","Repressed enhancer"))
length(JaenesL)
#2725 l1-l3
#all larval 3903

####### ecdf enhancer distance to fountains
JaenesL$distanceToFount<-mcols(distanceToNearest(JaenesL,fountains,ignore.strand=T))$distance

options(tibble.width=Inf)
dd1<-data.frame(JaenesL) %>%
  dplyr::group_by(topState_L3_chromHMM) %>%
  dplyr::mutate(ecd=ecdf(distanceToFount)(distanceToFount))

dd1$type<-factor(dd1$topState_L3_chromHMM, levels=c("Active enhancer" ,"Repressed enhancer","H3K27me3 repressed"))
levels(dd1$type)<-c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer")


# get some cutoffs
dd2<-dd1 %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(ecd10kb=ecdf(distanceToFount)(c(10000)),
                   ecd50kb=ecdf(distanceToFount)(c(50000)),
                   bracket1=distanceToFount[which.min(abs(ecd-0.6))],
                   bracket2=distanceToFount[which.min(abs(ecd-0.63))],
                   count=n())
dd1$type_count<-dd1$type
levels(dd1$type_count)<-paste0(c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer"),
       paste0(" (n=",dd2$count,")"))


p3<-ggplot(dd1, aes(x=distanceToFount/1000,y=ecd,color=type_count)) +
  geom_line(linewidth=0.9)+
  scale_color_manual(values=c("red","blue","darkgreen"))+
  xlab("Distance to fountain tip (kb)")+ylab("Cumulative distribution of distance")+
  theme(legend.position=c(0.7,0.5)) +
  coord_cartesian(xlim=c(1,150000/1000)) +
  geom_hline(yintercept=1,colour="darkgrey") +
  theme(legend.title=element_blank())+
  ggtitle("L1-L3 enhancers (Jaenes et al. 2018)")



# add 10kb quantile
p3<-p3+ geom_segment(x=10,y=dd2$ecd10kb[dd2$type=="Active enhancer"],xend=10,
                     yend=-Inf,linetype="dashed",colour="grey") +
  geom_segment(x=-Inf,y=dd2$ecd10kb[dd2$type=="Active enhancer"],
               xend=10,yend=dd2$ecd10kb[dd2$type=="Active enhancer"],
               linetype="dashed",colour="grey") +
  annotate(geom="text",x=11,y=0,label="10kb",color="darkgrey",hjust=0,size=3)+
  annotate(geom="text",x=0,y=dd2$ecd10kb[dd2$type=="Active enhancer"]*1.02,
           label=paste0(round((dd2$ecd10kb[dd2$type=="Active enhancer"])*100,0),"%"),
           color="darkgrey",vjust=0,size=3)
# add pvalues
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="Repressed enhancer"],alternative="less")
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")
ks.test(dd1$ecd[dd1$type=="Repressed enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")

p3<-p3+ geom_bracket(xmin=as.numeric(dd2[1,"bracket1"]/1000),
                     xmax=as.numeric(dd2[2,"bracket1"]/1000),
                     y.position=0.6, label="****",inherit.aes=F,
                     bracket.nudge.y=0.02, label.size=3,
                     tip.length=0.02,vjust=0.75) +
  geom_bracket(xmin=as.numeric(dd2[2,"bracket2"]/1000),
               xmax=as.numeric(dd2[3,"bracket2"]/1000),
               y.position=0.63, label="ns",inherit.aes=F,
               bracket.nudge.y=0.02, label.size=3,
               tip.length=0.02,vjust=0.2)

p3

# Number of enhancers per bin -------
binSize=6000
fountains$fountVcont<-"fountain tip"
nonFount$fountVcont<-"control"

allRegions<-c(fountains,nonFount)
allRegions<-resize(allRegions,width=binSize,fix="center")

allRegions$ActiveCount<-countOverlaps(allRegions,JaenesL[JaenesL$topState_L3_chromHMM=="Active enhancer"],ignore.strand=T,minoverlap=1)
allRegions$RepressedCount<-countOverlaps(allRegions,JaenesL[JaenesL$topState_L3_chromHMM=="Repressed enhancer"],ignore.strand=T,minoverlap=1)
allRegions$H3K27meRepressedCount<-countOverlaps(allRegions,JaenesL[JaenesL$topState_L3_chromHMM=="H3K27me3 repressed"],ignore.strand=T,minoverlap=1)


df<-data.frame(allRegions)

df1<-pivot_longer(df,cols=c("ActiveCount", "RepressedCount", "H3K27meRepressedCount"),
                  names_to="type",values_to="count")

df1$type<-gsub("Count$","",df1$type)
df1$type<-factor(df1$type,levels=c("Active", "Repressed","H3K27meRepressed"))
levels(df1$type)<-c("Active", "Repressed","H3K27me3")

p4<-ggplot(df1,aes(x=count,fill=type)) +
  geom_bar(position=position_dodge(preserve = "single")) +
  facet_wrap(.~fountVcont) +
  scale_fill_manual(values=c("red","blue","darkgreen")) +
  xlab(paste0("Number of enhancers per ",binSize/1000,"kb bin"))+
  ylab(paste0("Number of ",binSize/1000,"kb bins")) +
  ggtitle("L1-L3 enhancers (Jaenes et al. 2018)") +
  theme(legend.title=element_blank(), legend.position=c(0.8,0.6),
        legend.key.size = unit(0.3, 'cm'))


#p<-ggpubr::ggarrange(p1,p3,
#                     ggpubr::ggarrange(p2,p4,ncol=2,nrow=1),nrow=3,ncol=1)

p<-cowplot::plot_grid(p1,p2,p3,p4,nrow=2,ncol=2)
ggsave(paste0(finalFigDir,"/supplFigg_fountainsVenhancers.pdf"), p, device="pdf",
       width=19,height=18, units="cm")







