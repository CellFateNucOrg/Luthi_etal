library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(dplyr)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(plyranges)
library(rtracklayer)
library(tidyr)
library(ggtext)



theme_set(
  theme_classic(base_size=9)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=9),
          plot.title=ggtext::element_markdown(size=9),
          axis.title=ggtext::element_markdown(size=9),
          strip.text=ggtext::element_markdown(size=9),
          legend.title=ggtext::element_markdown(size=9))
)


projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseeq_DGE")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}


source(paste0(projectDir,"/functions.R"))
source(paste0(projectDir,"/functions_fountainPlots.R"))

fountains<-import(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.bed"))

chrRegions<-readRDS(paste0(publicDataDir,"/chr_regions_Rockman2009.RDS"))


chrdf<-data.frame(name=seqnames(Celegans),end=seqlengths(Celegans)/1e6)
chrdf<-chrdf[-7,]

fountdf<-data.frame(name=seqnames(fountains),start=start(fountains)/1e6,
                    end=end(fountains)/1e6)

chrdf$fountains<-table(fountdf$name)


# ideogram distribution by chromosome -----
p1<-ggplot(chrdf,aes(x=name,y=end))+
  geom_bar(stat="identity",width=0.4,fill="white",colour="black")+
  ylab("Position in Mb") +
  theme(axis.title.x=element_blank()) +
  geom_text(data=chrdf, aes(x=name, y=end,label=fountains),
            color="#00008b",vjust=-0.2,size=3)+
  geom_rect(data = fountdf, aes(xmin = as.numeric(name) - 0.2,
                                xmax = as.numeric(name) + 0.2,
                                ymax = end, ymin = start),
            col="#00008b66",linewidth=0.1)
#p1


# fountains per region vs length -----

regionsGR<-GRanges(chrRegions)
chrRegions$fountain_count<-countOverlaps(regionsGR,fountains)
chrRegions$seqnames<-factor(chrRegions$seqnames)
chrRegions$name<-factor(chrRegions$name)

## label left and right arms
armLabels<-chrRegions[chrRegions$name=="arm",]
armLabels$armLabel<-"right"
armLabels$armLabel[armLabels$end<8e6]<-"left"
armLabels$fountainPerMb<-armLabels$fountain_count/(armLabels$width/1e6)
armLabels$position<-ifelse(armLabels$fountainPerMb>mean(armLabels$fountainPerMb),"above","below")

p2<-ggplot(chrRegions,aes(x=width/1e6,y=fountain_count)) +
  geom_smooth(method="lm",color="black")+ geom_point(aes(color=seqnames,shape=name),size=3) +
  xlab("Size in Mb") + ylab("Number of fountains") +
  theme(legend.title=element_blank())

p2 <- p2 + geom_text_repel(data=subset(armLabels,position=="above"), aes(label=armLabel,color=seqnames),
                  max.overlaps = Inf,min.segment.length = 0, box.padding = 0.3,
                  nudge_y=5,nudge_x=-0.5,seed=42,size=3)
p2 <- p2 + geom_text_repel(data=subset(armLabels,position=="below"), aes(label=armLabel,color=seqnames),
                       max.overlaps = Inf,min.segment.length = 0, box.padding = 0.3,
                       nudge_y=-5,nudge_x=0.5,seed=42,size=3)
#p2


# Count fountains per Mb by chromosomal region ------
regionSize=1e6
gr<-GRanges(chrRegions)
tls<-tileGenome(seqlengths(Celegans),tilewidth = regionSize, cut.last.tile.in.chrom = T)
ol<-findOverlaps(tls,gr, minoverlap=0.501*regionSize)
mcols(tls)$name<-NA
mcols(tls)$name[queryHits(ol)] <- as.character(mcols(gr)$name[subjectHits(ol)])
tls1<-subsetByOverlaps(tls,gr, minoverlap=0.501*regionSize)
tls1<-tls1[width(tls1) == regionSize]
table(tls1$name)
#arm center    tip
#84     91      9

# get fountain count per Mb by region
tls1$fountain_count<-countOverlaps(tls1,fountains)

fountPerMb<-data.frame(tls1)
fountPerMb$seqnames<-factor(fountPerMb$seqnames)
fountPerMb$name<-factor(fountPerMb$name)


# wilcoxon by chr -----
stat.test <- fountPerMb %>%
  wilcox_test(fountain_count ~ seqnames) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "seqnames")
stat.test<-stat.test %>% dplyr::filter(p.adj.signif!="ns")
p3<-ggplot(fountPerMb,aes(x=seqnames,y=fountain_count)) +
  geom_beeswarm(alpha=0.6, fill="grey60", size=1) +
  stat_pvalue_manual(stat.test, label = "p.adj = {p.adj}",size=3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean,
               width = 0.4, color="red",size=1) +
  ylab(label=paste0("Fountains per Mb")) +
  theme(axis.title.x=element_blank())
#p3


# wilcoxon arms v center ------
stat.test <- fountPerMb %>%
  wilcox_test(fountain_count ~ name) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "name")
p4<-ggplot(fountPerMb,aes(x=name,y=fountain_count)) +
  geom_beeswarm(alpha=0.6, fill="grey60", size=1) +
  stat_pvalue_manual(stat.test, label = "{p.adj.signif}",size=3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean,
               width = 0.4, color="red",size=0.7) +
  theme(axis.title.x=element_blank()) +
  ylab(label=paste0("Fountains per Mb"))
#p4




#########################-
# rex sites and Albritton et al. loop anchors
#########################-


# enhancers
#rex<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))
rex<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/publicData/rexsites_Albritton2017_ce11.bed")
seqlevels(rex)<-seqlevels(Celegans)
seqinfo(rex)<-seqinfo(Celegans)
rexGap<-gaps(sort(resize(rex,width=1,fix="center")),ignore.strand=T)
rexGap<-rexGap[seqnames(rexGap)=="chrX"]
rexGap

fountX<-fountains[seqnames(fountains)=="chrX"]
colnames(mcols(fountX))[colnames(mcols(fountX))=="name"]<-"fountainName"

seqlevels(fountX)<-seqlevels(Celegans)
seqinfo(fountX)<-seqinfo(Celegans)

getRelativePosition<-function(rexGap,fountX){
  ol<-findOverlaps(rexGap,fountX)
  rexVfount=NULL
  for(i in unique(queryHits(ol))){
    fountInGap<-fountX[subjectHits(ol)[queryHits(ol)==i]]
    start<-start(rexGap[i])
    end<-end(rexGap[i])
    span<-end-start
    name<-paste0("rexGap",i)
    relPos<-round((start(fountInGap)-start+1000)/span,2)
    tmp<-data.frame(gapName=name,relPosition=relPos,fountainName=fountInGap$fountainName)
    if(is.null(rexVfount)){
      rexVfount<-tmp
    } else {
      rexVfount<-rbind(rexVfount,tmp)
    }
  }
  return(rexVfount)
}

rexVfount<-getRelativePosition(rexGap,fountX)
p5<-ggplot(rexVfount,aes(x=relPosition)) + geom_histogram(bins=100) +
  xlab("Relative position between adjacent rex sites") +
  ylab("Count")+
  geom_vline(xintercept=0,colour="red") +
  annotate("text",label="rex",x=0,y=6,color="red",hjust=-0.1)+
  geom_vline(xintercept=1,colour="red") +
  annotate("text",label="rex",x=1,y=6,color="red",hjust=1.1) +
  ggtitle("Distribution of fountains within xTADs") + geom_density(color="green")
#p5

####### ecdf enhancer distance to fountains
rex$distanceToFount<-mcols(distanceToNearest(rex,fountX,ignore.strand=T))$distance
rex$replicate="real"

rexRand<-rex
for(i in 1:20){
  newStarts<-sample(1:(17718942-2000), length(fountX), replace = F)
  fakeFount<-GRanges(seqnames="chrX",ranges=IRanges(start=newStarts,width=2000))
  seqlevels(fakeFount)<-seqlevels(Celegans)
  tmp<-rex
  tmp$replicate<-paste0("rand",i)
  tmp$distanceToFount<-mcols(distanceToNearest(rex,fakeFount,ignore.strand=T))$distance
  rexRand<-c(rexRand,tmp)
}

rexRand$type<-ifelse(rexRand$replicate=="real","real","rand")

options(tibble.width=Inf)
dd1<-data.frame(rexRand) %>%
  dplyr::group_by(replicate) %>%
  dplyr::mutate(ecd=ecdf(distanceToFount)(distanceToFount))
dd1$replicate


p6<-ggplot(dd1, aes(x=distanceToFount/1000,y=ecd,group_by=replicate,color=type)) +
  geom_line(linewidth=0.9)+
  xlab("Distance to nearest fountain tip (kb)")+ylab("Fraction")+
  scale_color_manual(labels=c("random fountains","real fountains"),values=c("grey","blue"))+
  coord_cartesian(xlim=c(1,150000/1000)) +
  geom_hline(yintercept=1,colour="darkgrey") +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.5)) +
  ggtitle("Random fountain position chosen from whole chromosome")
#p6


# choosing fountain starts for gene starts
genes<-readRDS(paste0(projectDir,"/ce11GeneGR_WS285.rds"))
genes<-genes[seqnames(genes)=="chrX"]

rexRand<-rex
for(i in 1:20){
  newStarts<-sample(start(genes[width(genes)>100]), length(fountX), replace = F)
  fakeFount<-GRanges(seqnames="chrX",ranges=IRanges(start=newStarts,width=2000))
  seqlevels(fakeFount)<-seqlevels(Celegans)
  tmp<-rex
  tmp$replicate<-paste0("rand",i)
  tmp$distanceToFount<-mcols(distanceToNearest(rex,fakeFount,ignore.strand=T))$distance
  rexRand<-c(rexRand,tmp)
}

rexRand$type<-ifelse(rexRand$replicate=="real","real","rand")

options(tibble.width=Inf)
dd1<-data.frame(rexRand) %>%
  dplyr::group_by(replicate) %>%
  dplyr::mutate(ecd=ecdf(distanceToFount)(distanceToFount))
dd1$replicate


p7<-ggplot(dd1, aes(x=distanceToFount/1000,y=ecd,group_by=replicate,color=type)) +
  geom_line(linewidth=0.9)+
  xlab("Distance to nearest fountain tip (kb)")+ylab("Fraction")+
  scale_color_manual(labels=c("random fountains","real fountains"),values=c("grey","blue"))+
  coord_cartesian(xlim=c(1,150000/1000)) +
  geom_hline(yintercept=1,colour="darkgrey") +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.5)) +
  ggtitle("Random fountain position chosen from starts of genes>100bp long")
#p7












# Assemble final figure ------
pnull<-NULL
pa<-cowplot::plot_grid(p1,p2,p4,p3,nrow=2,ncol=2,
                      labels=c("a ","b ","c ","d "),
                      align="h",rel_widths=c(0.7,1))
pb<-cowplot::plot_grid(p5,p6,pnull,p7,nrow=2,ncol=2,
                   labels=c("e ","f "," ","g "),
                   rel_widths=c(0.7,1))

p<-cowplot::plot_grid(pa,pb,nrow=2,ncol=1,align="v")

p<-annotate_figure(p, top = text_grob("Isiaka et al., Figure S2", size = 12))
ggsave(paste0(finalFigDir,"/FigS2_FountainsPer",regionSize/1e6,"Mb_armsVCenter.pdf"),p,width=18,height=27,
device="pdf",units="cm")

