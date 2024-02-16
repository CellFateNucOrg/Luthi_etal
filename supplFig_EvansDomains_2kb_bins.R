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
# fountains$AsymmetryTercile<-paste0("Asym.Q",fountains$AsymmetryTercile)
# fountains$AsymmetryQuintile<-paste0("Asym.Q",fountains$AsymmetryQuintile)





######################################
## chromatin domains ---------
######################################

domains<-import.bed(paste0(publicDataDir,"/chromDomains_L3_Evans2016_ce11.bed"))
seqlevels(domains)<-seqlevels(Celegans)
domains$name<-factor(domains$name)
domains$score<-factor(as.numeric(domains$name))

domainClrs<-c("#cc0000","#5b5b5b","#bcbcbc")

getDomainOLtable<-function(gr,domains){
  ol<-data.frame(findOverlaps(gr,domains,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$domain<-factor(domains$score[ol$subjectHits],levels=1:3)
  ol$width<-width(domains)[ol$subjectHits]
  df<-ol%>%group_by(domain) %>% summarise(domainFrequency=n(),domainWidth=sum(width))
  alldomains<-data.frame(domain=factor(1:3, levels=1:3))
  df<-left_join(alldomains,df)
  df[is.na(df)]<-0
  return(df)
}


tiles<-tileGenome(seqlengths(Celegans)[1:6],tilewidth=2000,cut.last.tile.in.chrom = T)
tiles<-getUpDownstreamFountainData(tiles,fountains)

tiles<-makeGRangesFromDataFrame(tiles,keep.extra.columns = T)

tiles$binnedDistance<-binByDistance1(tiles$distanceToFountain,maxDist=20000,binSize =2000,fullTile=T)
tiles<-flipDownstreamOrientation(tiles)
table(tiles$binnedDistance1,tiles$fountainLocation)
#tiles$AsymmetryTercile<-factor(tiles$AsymmetryTercile,levels=paste0("Asym.Q",1:3))
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
  group_by(binnedDistance1) %>%
  mutate(domainWidth=sum(domainWidth),domainFrequency=sum(domainFrequency)) %>%
  filter(binnedDistance1=="0") %>% select(-domain,-fountainLocation) %>% distinct()

p1<-ggplot(df,aes(x=binnedDistance1,y=domainWidth/1e7,fill=domain)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_manual(values=domainClrs) +
  xlab("Binned distance (kb)")+
  theme(legend.position = "none") +
  ggtitle(paste0("Chromatin domain width around autosomal fountain tips")) +
  ylab("Bin coverage x10 Mb")
p1<-p1 + geom_bar(data=boxdf,aes(x=binnedDistance1, y=domainWidth/1e7,fill=NULL),stat="identity",
                  alpha=0,color="black",linewidth=0.5)

p1


# genome average
domains$width=width(domains)
genomeWide<-data.frame(domains) %>% group_by(name, score) %>% summarise(domainWidth=sum(width))

p2<-ggplot(genomeWide,aes(x=1,y=domainWidth/1e7,fill=name)) +
  geom_bar(position="stack",stat="identity")  +
  scale_fill_manual(values=domainClrs) +
  xlab("x") + ylab("Genome-wide coverage x10 Mb") + ggtitle("")
p2

p<-ggarrange(p1,p2,nrow=1,ncol=2,widths=c(0.8,0.25))

p<-annotate_figure(p, top = text_grob("Isiaka et al., Supl Figure ", size = 14))
ggsave(paste0(finalFigDir,"/supplFig_EvansDomains_2kb_bins.pdf"), p, device="pdf",
       width=19,height=10, units="cm")

