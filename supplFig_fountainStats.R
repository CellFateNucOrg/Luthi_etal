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
library(GGally)

theme_set(
  theme_bw(base_size=9)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=9),
          plot.title=element_text(size=9),
          axis.title.y=element_text(size=9),
          axis.title.x=element_markdown(size=9),
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


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125c.rds"))
seqlevels(fountains)<-seqlevels(Celegans)
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
#rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")



# enhancers
#daughertyL3<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))

#metadata<-readRDS("/Users/semple/Documents/MeisterLab/papers/Isiaka_etal/wbGeneGR_WS285.rds")
#seqlevelsStyle(metadata)<-"UCSC"

fount<-data.frame(fountains)
borders<-quantile(fount$Asymetry,seq(0,1,0.2))
midpoints<-(borders[1:5] + borders[2:6])/2
p1<-ggplot(fount,aes(x=Asymetry)) +
  geom_histogram(bins=50,fill="grey80",color="black") +
  geom_vline(xintercept=borders,color="purple",linewidth=0.2) +
  xlab("Asymmetry") + ylab("Count") +
  annotate(geom="text",x=midpoints[1],y=160,label="bottom quintile\n(asymmetric left)",color="purple",size=3,lineheight = .8)+
  annotate(geom="text",x=midpoints[5],y=160,label="top quintile\n(asymmetric right)",color="purple",size=3,lineheight = .8)+
  annotate(geom="text",x=0.06,y=200,label="middle quintile\n(symmetric)",hjust=0,color="purple",size=3,lineheight = .8)+
  annotate(geom = "curve", x = 0.2, y = 215, xend = midpoints[3], yend = 200,
    curvature = .8, arrow = arrow(length = unit(1, "mm")))#+
  #ggtitle("Asymmetry score")
p1



#borders<-quantile(fount$symm._prom._fountain_score,seq(0,1,0.2))
#midpoints<-(borders[1:5] + borders[2:6])/2
med<-median(fount$symm._prom._fountain_score*1e3)
p2<-ggplot(fount,aes(x=symm._prom._fountain_score*1e3)) +
  geom_histogram(bins=50,fill="grey80",color="black") +
  geom_vline(xintercept=med,color="purple") +
  annotate(geom="text",x=med, y=95,
           label=paste0("median = ",round(med,2)),
           color="purple",size=3,hjust=-0.1)+
  xlab("Fountain score (x10<sup>-3</sup>)") + ylab("Count")
  #geom_vline(xintercept=borders,color="purple") +
  #annotate(geom="text",x=midpoints[1],y=85,label="Q1",color="purple",size=2.5)+
  #annotate(geom="text",x=midpoints[2],y=85,label="Q2",color="purple",size=2.5)+
  #annotate(geom="text",x=midpoints[3],y=85,label="Q3",color="purple",size=2.5)+
  #annotate(geom="text",x=midpoints[4],y=85,label="Q4",color="purple",size=2.5)+
  #annotate(geom="text",x=midpoints[5],y=85,label="Q5",color="purple",size=2.5)#+
  #ggtitle("Fountain score")
p2



#borders<-quantile(fount$length_bp/1000,seq(0,1,0.2))
#midpoints<-(borders[1:5] + borders[2:6])/2
med<-median(fount$length_bp/1e3)
p4<-ggplot(fount,aes(x=length_bp/1e3)) +
  geom_histogram(bins=50,fill="grey80",color="black") +
  xlab("Fountain length (kb)") + ylab("Count") +
  geom_vline(xintercept=med,color="purple") +
  annotate(geom="text",x=med,y=85,
           label=paste0("median = ",round(med,2)),
           color="purple",size=3,hjust=-0.1)
  # geom_vline(xintercept=borders,color="purple") +
  # annotate(geom="text",x=midpoints[1],y=200,label="Q1",color="purple",size=2.5)+
  # annotate(geom="text",x=midpoints[2],y=200,label="Q2",color="purple",size=2.5)+
  # annotate(geom="text",x=midpoints[3],y=200,label="Q3",color="purple",size=2.5)+
  # annotate(geom="text",x=midpoints[4],y=200,label="Q4",color="purple",size=2.5)+
  # annotate(geom="text",x=midpoints[5],y=200,label="Q5",color="purple",size=2.5)#+
  #ggtitle("Fountain length")
p4


p6<-ggplot(fount,aes(x=symm._prom._fountain_score*1e3,y=coh1)) + geom_point() +
  geom_smooth(method="lm") +
  coord_cartesian(ylim=c(min(fount$coh1),max(fount$coh1)))+
  stat_cor(method = "pearson", label.x = min(fount$symm._prom._fountain_score*1e3),
           label.y = max(fount$coh1)*0.95,size=3)+
  ylab("COH-1 ChIP at tip") + xlab("Fountain score (x10<sup>-3</sup>)")
p6






p<-cowplot::plot_grid(p1,p2,p4,p6,nrow=2,ncol=2,align="v",
                      labels=c("g ","h ","i ", "j "))
p<-annotate_figure(p, top = text_grob("Lüthi et al., Figure S1", size = 12))
ggsave(paste0(finalFigDir,"/supplFig_fountainStats.pdf"), p, device="pdf",
       width=18,height=17, units="cm")



#export(fountains[fountains$AsymmetryQuintile==1],"fountains_assymmetryQ1.bed")
#export(fountains[fountains$AsymmetryQuintile==3],"fountains_assymmetryQ3.bed")
#export(fountains[fountains$AsymmetryQuintile==5],"fountains_assymmetryQ5.bed")
