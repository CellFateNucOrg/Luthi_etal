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
source(paste0(projectDir,"/functions_plotting.R"))


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125c.rds"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")


# control bins
nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=2000,fix="center")
rtracklayer::export(nonFount,paste0(fountainsDir,"/controlsForFountains_base0_uncorrected_20240125.bed"))

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
levels(dd1$type)<-c("Active enhancer", "Repressed enhancer", "H3K27me3 region")

table(dd1$type)
# to get some stats
dd2<-dd1 %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(ecd10kb=ecdf(distanceToFount)(c(10000)),
                   ecd50kb=ecdf(distanceToFount)(c(50000)),
                   bracket1=distanceToFount[which.min(abs(ecd-0.6))],
                   bracket2=distanceToFount[which.min(abs(ecd-0.63))],
                   count=n())

dd1$type_count<-dd1$type
levels(dd1$type_count)<-paste0(c("Active enhancer", "Repressed enhancer", "H3K27me3 region"),
                               paste0(" (n=",dd2$count,")"))

p1<-ggplot(dd1, aes(x=distanceToFount/1000,y=ecd,color=type_count)) +
  geom_line(linewidth=0.9)+
  scale_color_manual(values=c("red","blue","darkgreen"))+
  xlab("Distance to fountain 2kb tip bin (kb)")+ylab("Cumulative distribution of distance")+
  coord_cartesian(xlim=c(1,150000/1000)) +
  geom_hline(yintercept=1,colour="darkgrey") +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.5)) +
  ggtitle("L3 enhancers (Daugherty et al. 2017)")


# add 10kb quantile
p1<-p1+ geom_segment(x=10,y=dd2$ecd10kb[dd2$type=="Active enhancer"],xend=10,
                     yend=-Inf,linetype="dashed",colour="grey") +
  geom_segment(x=-Inf,y=dd2$ecd10kb[dd2$type=="Active enhancer"],
               xend=10,yend=dd2$ecd10kb[dd2$type=="Active enhancer"],
               linetype="dashed",colour="grey") +
  annotate(geom="text",x=11,y=0,label="10kb",color="darkgrey",hjust=0,size=3)+
  annotate(geom="text",x=1,y=dd2$ecd10kb[dd2$type=="Active enhancer"]*1.02,
           label=paste0(round((dd2$ecd10kb[dd2$type=="Active enhancer"])*100,0),"%"),
           color="darkgrey",vjust=0,size=3)


pActive<-ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="Repressed enhancer"],alternative="less")$p.value
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="H3K27me3 region"],alternative="less")
pRep<-ks.test(dd1$ecd[dd1$type=="Repressed enhancer"],dd1$ecd[dd1$type=="H3K27me3 region"],alternative="less")$p.value
pActive
#pActive<-formatCustomSci(pActive)
pActive<-ifelse(pActive<0.0001,"****","***")
pRep<-ifelse(pRep>0.05,"ns",pRep)


p1<-p1+ geom_bracket(xmin=as.numeric(dd2[1,"bracket1"]/1000),
                     xmax=as.numeric(dd2[2,"bracket1"]/1000),
                     y.position=0.6, label=pActive,inherit.aes=F, bracket.nudge.y=0.02,
                     tip.length=0.02,vjust=0.75,hjust=0.4,label.size=3) +
  geom_bracket(xmin=as.numeric(dd2[2,"bracket2"]/1000),
               xmax=as.numeric(dd2[3,"bracket2"]/1000),
               y.position=0.63, label=pRep,inherit.aes=F, bracket.nudge.y=0.02,
               tip.length=0.02,vjust=0.2,hjust=0.4,label.size=3)

p1


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



df1$fountVcont<-factor(df1$fountVcont, levels=c("control","fountain tip"))
levels(df1$fountVcont)<-c("control","fountain\ntip")

p2<-ggplot(df1,aes(x=fountVcont,group=count,fill=factor(count))) +
  geom_bar(position = position_fill(reverse = TRUE),color="black",linewidth=0.1) +
  facet_wrap(.~type) +
  scale_fill_grey(paste0("Enhancers\nper ",binSize/1000,"kb bin"),start=1,end=0) +
  xlab(paste0(""))+
  ylab(paste0("Fraction of ",binSize/1000,"kb bins"))+
  theme(legend.title=element_text(size=8,angle=0),
        legend.position="right",
        legend.key.size = unit(0.3, 'cm'))+
  #guides(fill=guide_legend(title.position="left",title.vjust=0))+
  ggtitle("L3 enhancers (Daugherty et al. 2017)")
p2


binSize=6000
allRegions<-c(fountains)
allRegions<-resize(allRegions,width=binSize,fix="center")

allRegions$ActiveCount<-countOverlaps(allRegions,daughertyL3[daughertyL3$L3_chromHMMState=="L3_activeEnhancer"],ignore.strand=T,minoverlap=1)


tiles<-tileGenome(seqlengths=seqlengths(Celegans),tilewidth=binSize,cut.last.tile.in.chrom = T)
tiles$ActiveCount<-countOverlaps(tiles,daughertyL3[daughertyL3$L3_chromHMMState=="L3_activeEnhancer"],ignore.strand=T,minoverlap=1)

tileSummary<-data.frame(tiles) %>%
  dplyr::group_by(ActiveCount) %>% summarise(bins=n())

df2<-data.frame(allRegions) %>%
  dplyr::group_by(ActiveCount) %>% summarise(bins=n())

df3<-left_join(tileSummary,df2,by=join_by(ActiveCount),suffix=c(".all",".fount"))
df3[is.na(df3)]<-0
df3$bins.nonFount<-df3$bins.all-df3$bins.fount


df4<-pivot_longer(df3,cols=c(bins.fount,bins.nonFount),names_to="fountVnonFount",
                  values_to="count")


p3<-ggplot(df4,aes(x=factor(ActiveCount),y=count,fill=fountVnonFount)) +
  geom_bar(stat = "identity",position = position_fill(reverse=T),
           color="black",width=0.8)+
  scale_fill_grey(paste0(binSize/1000,"kb bin type"),start=0,end=1,
                  labels=c("Fountain","Not fountain")) +
  geom_hline(yintercept=length(fountains)/length(tiles),color="red",
             linetype="dashed") +
  geom_text(aes(x=factor(ActiveCount),y=1.05,label=bins.all),color="blue",
            size=2.5) +
  xlab(paste0("Number of active enhancers per ",binSize/1000,"kb bin")) +
  ylab(paste0("Fraction of ",binSize/1000,"kb bins")) +
  theme(legend.position="bottom",legend.key.size = unit(0.3, 'cm'))+
  ggtitle("L3 enhancers (Daugherty et al. 2017)")


df<-data.frame(allRegions)
df$ActiveCount<-factor(df$ActiveCount)
stat.test<-df  %>% wilcox_test(symm._prom._fountain_score~ActiveCount,ref.group="0") %>%
  add_xy_position() %>% #p_format(new.col=T,accuracy=1e-32)%>%
  mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")


df1<-df%>% group_by(ActiveCount) %>% summarise(total=n())
df1

p4<-ggplot(df,aes(x=factor(ActiveCount),y=symm._prom._fountain_score*1e3)) +
  geom_boxplot(notch=T,outlier.shape=NA,fill="grey80") +
  xlab("Number of active enhancers per 6kb bin") +
  ylab("Fountain prominence score (x10<sup>-3</sup>)")+
  #stat_pvalue_manual(stat.test, label = "p.adj.format",remove.bracket=T,
  #                   y=2.3) +
  geom_signif(y_position=2.3,
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0)+
  geom_text(aes(x=ActiveCount,y=0,label=total),data=df1,color="blue",
            size=2.5)+
  ggtitle("L3 enhancers (Daugherty et al. 2017)")
p4






#########################-
# Jaenes ATAC enhancers -----
#########################-


# enhancers
JaenesL<-readRDS(paste0(publicDataDir,"/Jaenes2018_enhancers_ce11_stages_chromHMM.rds"))
JaenesL<-JaenesL%>% filter(topState_L3_chromHMM %in% c("Active enhancer","H3K27me3 repressed","Repressed enhancer"))
#7543
#2725 l1-l3
#all larval 3903

####### ecdf enhancer distance to fountains
JaenesL$distanceToFount<-mcols(distanceToNearest(JaenesL,fountains,ignore.strand=T))$distance

options(tibble.width=Inf)
dd1<-data.frame(JaenesL) %>%
  dplyr::group_by(topState_L3_chromHMM) %>%
  dplyr::mutate(ecd=ecdf(distanceToFount)(distanceToFount))

dd1$type<-factor(dd1$topState_L3_chromHMM, levels=c("Active enhancer" ,"Repressed enhancer","H3K27me3 repressed"))
levels(dd1$type)<-c("Active enhancer", "Repressed enhancer", "H3K27me3 region")

table(dd1$type)

# get some cutoffs
dd2<-dd1 %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(ecd10kb=ecdf(distanceToFount)(c(10000)),
                   ecd50kb=ecdf(distanceToFount)(c(50000)),
                   bracket1=distanceToFount[which.min(abs(ecd-0.6))],
                   bracket2=distanceToFount[which.min(abs(ecd-0.63))],
                   count=n())

dd1$type_count<-dd1$type
levels(dd1$type_count)<-paste0(c("Active enhancer", "Repressed enhancer", "H3K27me3 region"),
       paste0(" (n=",dd2$count,")"))


p5<-ggplot(dd1, aes(x=distanceToFount/1000,y=ecd,color=type_count)) +
  geom_line(linewidth=0.9)+
  scale_color_manual(values=c("red","blue","darkgreen"))+
  xlab("Distance to fountain 2kb tip bin (kb)")+ylab("Cumulative distribution of distance")+
  theme(legend.position=c(0.7,0.5)) +
  coord_cartesian(xlim=c(1,150000/1000)) +
  geom_hline(yintercept=1,colour="darkgrey") +
  theme(legend.title=element_blank())+
  ggtitle("L3 enhancers (Jaenes et al. 2018)")



# add 10kb quantile
p5<-p5+ geom_segment(x=10,y=dd2$ecd10kb[dd2$type=="Active enhancer"],xend=10,
                     yend=-Inf,linetype="dashed",colour="grey") +
  geom_segment(x=-Inf,y=dd2$ecd10kb[dd2$type=="Active enhancer"],
               xend=10,yend=dd2$ecd10kb[dd2$type=="Active enhancer"],
               linetype="dashed",colour="grey") +
  annotate(geom="text",x=11,y=0,label="10kb",color="darkgrey",hjust=0,size=3)+
  annotate(geom="text",x=1,y=dd2$ecd10kb[dd2$type=="Active enhancer"]*1.02,
           label=paste0(round((dd2$ecd10kb[dd2$type=="Active enhancer"])*100,0),"%"),
           color="darkgrey",vjust=0,size=3)


pActive<-ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="Repressed enhancer"],alternative="less")$p.value
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="H3K27me3 region"],alternative="less")
pRep<-ks.test(dd1$ecd[dd1$type=="Repressed enhancer"],dd1$ecd[dd1$type=="H3K27me3 region"],alternative="less")$p.value

#pActive<-formatCustomSci(pActive)
pActive<-ifelse(pActive<0.0001,"****","***")
pRep<-ifelse(pRep>0.05,"ns",pRep)


p5<-p5+ geom_bracket(xmin=as.numeric(dd2[1,"bracket1"]/1000),
                     xmax=as.numeric(dd2[2,"bracket1"]/1000),
                     y.position=0.6, label=pActive, inherit.aes=F,
                     bracket.nudge.y=0.02, label.size=3,
                     tip.length=0.02,vjust=0.75,hjust=0.4) +
  geom_bracket(xmin=as.numeric(dd2[2,"bracket2"]/1000),
               xmax=as.numeric(dd2[3,"bracket2"]/1000),
               y.position=0.63, label=pRep,inherit.aes=F,
               bracket.nudge.y=0.02, label.size=3,
               tip.length=0.02,vjust=0.2,hjust=0.4)

p5

# Number of enhancers per bin -------
binSize=6000
fountains$fountVcont<-"fountain tip"
nonFount$fountVcont<-"control"
tiles<-tileGenome(seqlengths=seqlengths(Celegans),tilewidth=binSize,cut.last.tile.in.chrom = T)

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

df1$fountVcont<-factor(df1$fountVcont, levels=c("control","fountain tip"))
levels(df1$fountVcont)<-c("control","fountain\ntip")



p6<-ggplot(df1,aes(x=fountVcont,group=count,fill=factor(count))) +
  geom_bar(position = position_fill(reverse = TRUE),color="black",linewidth=0.1) +
  facet_wrap(.~type) +
  scale_fill_grey(paste0("Enhancers\nper ",binSize/1000,"kb bin"),start=1,end=0) +
  xlab(paste0(""))+
  ylab(paste0("Fraction of ",binSize/1000,"kb bins"))+
  theme(legend.title=element_text(size=8,angle=0),
        legend.position="right",
        legend.key.size = unit(0.3, 'cm'))+
  ggtitle("L3 enhancers (Jaenes et al. 2018)")
p6


getRandomClustering<-function(allRegions,numSim=100){
  fountRegions<-allRegions[allRegions$fountVcont=="fountain tip"]
  numEnh<-sum(fountRegions$ActiveCount)
  numRegions<-length(fountRegions)
  mat<-matrix(data=0,nrow=numSim,ncol=30)
  colnames(mat)<-0:29
  for(sim in 1:numSim){
    s<-sample(1:numRegions,numEnh,replace=T)
    d<-table(table(s))
    mat[sim,"0"]<-numRegions-sum(d)
    mat[sim,names(d)]<-d
  }
  return(mat)
}


binSize=6000
allRegions<-c(fountains)
allRegions<-resize(allRegions,width=binSize,fix="center")

allRegions$ActiveCount<-countOverlaps(allRegions,JaenesL[JaenesL$topState_L3_chromHMM=="Active enhancer"],ignore.strand=T,minoverlap=1)


tiles<-tileGenome(seqlengths=seqlengths(Celegans),tilewidth=binSize,cut.last.tile.in.chrom = T)
tiles$ActiveCount<-countOverlaps(tiles,JaenesL[JaenesL$topState_L3_chromHMM=="Active enhancer"],ignore.strand=T,minoverlap=1)

tileSummary<-data.frame(tiles) %>%
  dplyr::group_by(ActiveCount) %>% summarise(bins=n())

df2<-data.frame(allRegions) %>%
  dplyr::group_by(ActiveCount) %>% summarise(bins=n())

df3<-left_join(tileSummary,df2,by=join_by(ActiveCount),suffix=c(".all",".fount"))

df3$bins.nonFount<-df3$bins.all-df3$bins.fount

df4<-pivot_longer(df3,cols=c(bins.fount,bins.nonFount),names_to="fountVnonFount",
             values_to="count")


p7<-ggplot(df4,aes(x=factor(ActiveCount),y=count,fill=fountVnonFount)) +
  geom_bar(stat = "identity",position = position_fill(reverse=T),
           color="black",width=0.8)+
  scale_fill_grey(paste0(binSize/1000,"kb bin type"),start=0,end=1,
                  labels=c("Fountain","Not fountain")) +
  geom_hline(yintercept=length(fountains)/length(tiles),color="red",
             linetype="dashed") +
  geom_text(aes(x=factor(ActiveCount),y=1.05,label=bins.all),color="blue",
            size=2.5) +
  xlab(paste0("Number of active enhancers per ",binSize/1000,"kb bin")) +
  ylab(paste0("Fraction of ",binSize/1000,"kb bins")) +
  theme(legend.position="bottom",legend.key.size = unit(0.3, 'cm'))+
  ggtitle("L3 enhancers (Jaenes et al. 2018)")


df<-data.frame(allRegions)
df$ActiveCount<-factor(df$ActiveCount)
stat.test<-df  %>% wilcox_test(symm._prom._fountain_score~ActiveCount,ref.group="0") %>%
  add_xy_position() %>%
  mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test
df1<-df%>% group_by(ActiveCount) %>% summarise(total=n())
df1

p8<-ggplot(df,aes(x=factor(ActiveCount),y=symm._prom._fountain_score*1e3)) +
  geom_boxplot(notch=T,outlier.shape=NA,fill="grey80") +
  xlab("Number of active enhancers per 6kb bin") +
  ylab("Fountain prominence score (x10<sup>-3</sup>)")+
  geom_signif(y_position=rep(c(2.2,2.35), nrow(stat.test))[1:nrow(stat.test)],
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0,step_increase=0.1)+
  geom_text(aes(x=ActiveCount,y=0,label=total),data=df1,color="blue",size=2.5)+
  ggtitle("L3 enhancers (Jaenes et al. 2018)")
p8




p<-cowplot::plot_grid(p1,p5,p2,p6,p3,p7,p4,p8,nrow=4, ncol=2,
                      labels=c("a ","b ","c ", "d ","e ","f ","g ","h "),
                      align = "v", axis="tb")
p<-annotate_figure(p, top = text_grob("Lüthi et al., Figure S3", size = 12))
ggsave(paste0(finalFigDir,"/FigS3_fountainsVenhancers.pdf"), p, device="pdf",
       width=19,height=29, units="cm")


