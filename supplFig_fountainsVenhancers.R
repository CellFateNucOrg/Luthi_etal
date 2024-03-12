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


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
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
levels(dd1$type)<-c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer")

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
levels(dd1$type_count)<-paste0(c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer"),
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
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")
pRep<-ks.test(dd1$ecd[dd1$type=="Repressed enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")$p.value
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

# p2<-ggplot(df1,aes(x=count,fill=type)) +
#   geom_bar(position=position_dodge(preserve = "single")) +
#   facet_wrap(.~fountVcont) +
#   scale_fill_manual(values=c("red","blue","darkgreen")) +
#   xlab(paste0("Number of enhancers per ",binSize/1000,"kb bin"))+
#   ylab(paste0("Number of ",binSize/1000,"kb bins"))+
#   theme(legend.title=element_blank(), legend.position=c(0.8,0.7),
#         legend.key.size = unit(0.3, 'cm'))+
#   ggtitle("L3 enhancers (Daugherty et al. 2017)")

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

mat<-getRandomClustering(allRegions,numSim=1000)
randStats<-rbind(colMeans(mat),apply(mat,2,quantile,probs=c(0.025, 0.975)))
randStats<-randStats[,1:(max(df1$count)+1)]
rownames(randStats)<-c("mean","lci","uci")
randdf<-data.frame(t(randStats))
randdf$count<-rownames(randdf)
randdf$data<-"expected"


df2<-df1 %>% dplyr::filter(fountVcont=="fountain\ntip",type=="Active") %>%
  dplyr::group_by(count) %>% summarise(mean=n(),lci=NA,uci=NA)
#df2$count<-factor(df2$count)
df2$data<-"observed"

randdf<-rbind(randdf,df2)
randdf$count<-factor(randdf$count,levels=0:max(df1$count))

p3<-ggplot(randdf,aes(x=factor(count),y=mean)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.3,linewidth=0.5) +
  geom_point(aes(color=data,shape=data),size=2) +
  scale_shape_manual(values=c(3,16)) +
  xlab("Number of active enhancers for 6kb fountain tip bin") +
  ylab("Number of 6kb fountain tip bins")+
  scale_color_manual(values=c("black","red"))+
  theme(legend.position=c(0.8,0.6),legend.title=element_blank())+
  ggtitle("L3 enhancers (Daugherty et al. 2017)")
p3

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
levels(dd1$type)<-c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer")

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
levels(dd1$type_count)<-paste0(c("Active enhancer", "Repressed enhancer", "H3K27me3 enhancer"),
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

# add pvalues
# formatCustomSci <- function(x) {     # Create user-defined function
#   x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
#   alpha <- round(as.numeric(x_sci[ , 1]),1)
#   power <- as.integer(x_sci[ , 2])
#   if(x!=0){
#     pval<-paste0(alpha,"x10",power)
#   } else {
#     pval<-"<2.2x10-16"
#   }
#   return(pval)
# }

pActive<-ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="Repressed enhancer"],alternative="less")$p.value
ks.test(dd1$ecd[dd1$type=="Active enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")
pRep<-ks.test(dd1$ecd[dd1$type=="Repressed enhancer"],dd1$ecd[dd1$type=="H3K27me3 enhancer"],alternative="less")$p.value

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

# p4<-ggplot(df1,aes(x=count,fill=type)) +
#   geom_bar(position=position_dodge(preserve = "single")) +
#   facet_wrap(.~fountVcont) +
#   scale_fill_manual(values=c("red","blue","darkgreen")) +
#   xlab(paste0("Number of enhancers per ",binSize/1000,"kb bin"))+
#   ylab(paste0("Number of ",binSize/1000,"kb bins")) +
#   ggtitle("L1-L3 enhancers (Jaenes et al. 2018)") +
#   theme(legend.title=element_blank(), legend.position=c(0.8,0.6),
#         legend.key.size = unit(0.3, 'cm'))

p6<-ggplot(df1,aes(x=fountVcont,group=count,fill=factor(count))) +
  geom_bar(position = position_fill(reverse = TRUE),color="black",linewidth=0.1) +
  facet_wrap(.~type) +
  scale_fill_grey(paste0("Enhancers\nper ",binSize/1000,"kb bin"),start=1,end=0) +
  xlab(paste0(""))+
  ylab(paste0("Fraction of ",binSize/1000,"kb bins"))+
  theme(legend.title=element_text(size=8,angle=0),
        legend.position="right",
        legend.key.size = unit(0.3, 'cm'))+
  #guides(fill=guide_legend(title.position="left",title.vjust=0,title.hjust=1))+
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

mat<-getRandomClustering(allRegions,numSim=1000)
randStats<-rbind(colMeans(mat),apply(mat,2,quantile,probs=c(0.025, 0.975)))
randStats<-randStats[,1:(max(df1$count)+1)]
rownames(randStats)<-c("mean","lci","uci")
randdf<-data.frame(t(randStats))
randdf$count<-rownames(randdf)
randdf$data<-"expected"


df2<-df1 %>% dplyr::filter(fountVcont=="fountain\ntip",type=="Active") %>%
  dplyr::group_by(count) %>% summarise(mean=n(),lci=NA,uci=NA)
#df2$count<-factor(df2$count,levels=as.character(0:max(df1$count)))
df2$data<-"observed"

randdf<-rbind(randdf,df2)
randdf$count<-factor(randdf$count,levels=as.character(0:max(df1$count)))

p7<-ggplot(randdf,aes(x=count,y=mean)) +
  geom_errorbar(aes(ymin=lci, ymax=uci),width=0.4,linewidth=0.5) +
  geom_point(aes(color=data,shape=data),size=2) +
  scale_shape_manual(values=c(3,16)) +
  xlab("Number of active enhancers for 6kb fountain tip bin") +
  ylab("Number of 6kb fountain tip bins")+
  scale_color_manual(values=c("black","red"))+
  theme(legend.position=c(0.8,0.6),legend.title=element_blank()) +
  ggtitle("L3 enhancers (Jaenes et al. 2018)")
p7




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

#ggplot(data=df[df$ActiveCount!=0,],aes(x=ActiveCount,y=log2FoldChange,fill=ActiveCount)) + geom_beeswarm() +
#  geom_hline(yintercept=0) + scale_fill_grey(start=1,end=0.4)


#' Function for adding count data to plots
count_data <- function (y,ymax=5){
  df <- data.frame(y = ymax, label = length(y))
  return(df)
}


p8<-ggplot(data=df,aes(x=ActiveCount,y=log2FoldChange)) +
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
p8





#p<-ggpubr::ggarrange(p1,p3,
#                     ggpubr::ggarrange(p2,p4,ncol=2,nrow=1),nrow=3,ncol=1)

p<-cowplot::plot_grid(p1,p5,p2,p6,p3,p7,p4,p8,nrow=4,ncol=2,
                      rel_widths=c(1,1,1,1,1,1,1,1),
                      labels=c("a ","b ","c ", "d ","e ","f ","g ","h "),
                      align = "v", axis="tb")
p<-annotate_figure(p, top = text_grob("Isiaka et al., Supl. Figure", size = 14))
ggsave(paste0(finalFigDir,"/supplFig_fountainsVenhancers.pdf"), p, device="pdf",
       width=19,height=29, units="cm")


