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


theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(),
          axis.title.x=ggtext::element_markdown())
)


projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseeq_DGE")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

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
  theme_classic() + ylab("Position in Mb") + xlab("") +
  geom_text(data=chrdf, aes(x=name, y=end,label=fountains),color="#00008b",vjust=-0.2)+
  geom_rect(data = fountdf, aes(xmin = as.numeric(name) - 0.2,
                                xmax = as.numeric(name) + 0.2,
                                ymax = end, ymin = start),
            col="#00008b66",linewidth=0.1)
p1


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

p4<-ggplot(chrRegions,aes(x=width/1e6,y=fountain_count)) +  theme_classic() +
  geom_smooth(method="lm",color="black")+ geom_point(aes(color=seqnames,shape=name),size=4) +
  xlab("Size in Mb") + ylab("Number of fountains")

p4 <- p4 + geom_text_repel(data=subset(armLabels,position=="above"), aes(label=armLabel,color=seqnames),
                  max.overlaps = Inf,min.segment.length = 0, box.padding = 0.3,
                  nudge_y=5,nudge_x=-0.5,seed=42)
p4 <- p4 + geom_text_repel(data=subset(armLabels,position=="below"), aes(label=armLabel,color=seqnames),
                       max.overlaps = Inf,min.segment.length = 0, box.padding = 0.3,
                       nudge_y=-5,nudge_x=0.5,seed=42)
p4


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
p2<-ggplot(fountPerMb,aes(x=seqnames,y=fountain_count)) +
  geom_beeswarm(alpha=0.6, fill="grey60", size=1) +
  theme_classic() +
  stat_pvalue_manual(stat.test, label = "p.adj = {p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.4,
               color="red",size=1) + xlab("") +
  ylab(label=paste0("Fountains per Mb")) #+ ggtitle("Fountains per Mb by chromosome")
p2


# wilcoxon arms v center ------
stat.test <- fountPerMb %>%
  wilcox_test(fountain_count ~ name) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "name")
p3<-ggplot(fountPerMb,aes(x=name,y=fountain_count)) +
  geom_beeswarm(alpha=0.6, fill="grey60", size=1) +
  theme_classic() +
  stat_pvalue_manual(stat.test, label = "{p.adj.signif}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .4,
               color="red",size=1) + xlab("") +
  ylab(label=paste0("Fountains per Mb"))
p3


# Assemble final figure ------
p<-ggarrange(p4,p1,p2,p3,nrow=2,ncol=2, labels=c("A ","C ","B ","D "),widths=c(1,0.6,1,0.6))
p<-annotate_figure(p, top = text_grob("Isiaka et al., Supplementary Figure ", size = 14))
ggsave(paste0(finalFigDir,"/FountainsPer",regionSize/1e6,"Mb_armsVCenter.pdf"),p,width=19,height=19,
device="pdf",units="cm")

