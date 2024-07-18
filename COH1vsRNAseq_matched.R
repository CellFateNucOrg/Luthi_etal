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

coh1<-import("/Volumes/MeisterLab/publicData/ChIPseq_realigned/modEncode_SMC/chip-seq-pipeline2/bigwig_fc/COH1_fem2_AD_SDQ0809_rep.pooled_x_ctl.pooled.fc.signal.bigwig")



########
tiles<-tileGenome(seqlengths(Celegans),tilewidth=2000,cut.last.tile.in.chrom = T)

fdist<-distanceToNearest(tiles,fountains,ignore.strand=T)
tiles$fountDistance_kb<-NA
tiles[queryHits(fdist)]$fountDistance_kb<-mcols(fdist)$distance/1000
tiles<-tiles[!is.na(tiles$fountDistance_kb)] # get rid of mtDNA

cov<-coverage(coh1,weight="score")

tiles<-binnedAverage(tiles,cov[seqlevels(Celegans)],varname="COH1",na.rm=T)

# COH-1 signal vs fountain distance
minDistance=30
minCOH1=1
df<-data.frame(tiles)


p1a<-ggplot(df, aes(x=fountDistance_kb,y=COH1)) +
  geom_bin2d(bins=100) +
  xlab("Distance to nearest fountain (kb)") +
  ylab("COH-1 ChIP signal per 2kb bin") +
  scale_fill_gradient2(low="lightgrey", mid = "darkblue",high="lightblue",midpoint=100)

p2<-p1a+annotate("rect", xmin = 0, xmax = 3 , ymin = minCOH1, ymax = max(nonfountbins$COH1),
           alpha = .1, colour = "red", fill=NA) +
  annotate("text", x=0, y=max(nonfountbins$COH1)*1.2, label="fountain\ntip bins\nsampled",
           hjust=0.2,vjust=0.3) +
  annotate("rect", xmin = minDistance, xmax = 200 , ymin = minCOH1, ymax = max(nonfountbins$COH1),
           alpha = .1, colour = "red",fill=NA)+
annotate("text", x=100, y=max(nonfountbins$COH1)*1.1, label=paste0("non fountain bins sampled"))
p2

nonfountbins<-df[df$COH1>minCOH1 & df$fountDistance_kb>minDistance,]
dim(nonfountbins) #1392

fountbins<-df[df$COH1>minCOH1 & df$COH1<=max(nonfountbins$COH1) & df$fountDistance_kb<3,]
dim(fountbins) #1915


set1=fountbins
set2=nonfountbins
numQuantiles=10
metricName="COH1"

#' Function to create two sets of elements with matching metric
#'
#' @param set1 dataframe of first set
#' @param set2 dataframe of first set
#' @param metricName name of column containing the metric to match between the
#' two sets
#' @param numQuantiles number of quantiles to use for matching the metric
#' @retrun dataframe with the selected matched elements of each set.
#' @export
getMatchingSets<-function(set1,set2,metricName,numQuantiles=10){
  set1$set<-1
  set2$set<-2
  all<-rbind(set1,set2)
  breaks<-quantile(all[,metricName],probs=seq(0,1,1/(numQuantiles)))
  all$quantile<-cut(all[,metricName],breaks=breaks,labels=1:numQuantiles,include.lowest=T)
  newall<-NULL
  for(q in 1:numQuantiles){
    tb<-table(all[all$quantile==q,"set"])
    #print(paste0("quantile ",q," has:"))
    #print(tb)
    smallset<-as.numeric(which.min(tb))
    smallsetsize<-as.numeric(min(tb))
    bigset<-as.numeric(which.max(tb))
    bigsetsize<-as.numeric(max(tb))
    # take all of small set
    if(is.null(newall)){
      newall<-all[all$quantile==q & all$set==smallset,]
    } else {
      newall<-rbind(newall,all[all$quantile==q & all$set==smallset,])
    }
    #sample same number from bigger set
    idx<-sample(1:bigsetsize,size=smallsetsize,replace=F)
    #print(paste(smallsetsize,bigsetsize,length(idx)))
    tmp<-all[all$quantile==q & all$set==bigset,]
    newall<-rbind(newall,tmp[idx,])
  }
return(newall)
}

seed=1243
set.seed(seed)
dd<-getMatchingSets(fountbins,nonfountbins,metricName="COH1")

table(dd$set)

dd$set<-factor(dd$set)
levels(dd$set)<-c("fountain tips +-2kb",paste0(">",minDistance,"kb from fountain tips"))

p3<-ggplot(dd,aes(x=COH1)) + geom_histogram() +
  facet_wrap(~set) +
  geom_text(data = dd %>% group_by(set) %>% summarise(count=n()),
            aes(label = paste("Count:",count), y = Inf, x  = Inf), vjust = 1.2, hjust = 1)+
  ylab("Number of 2kb bins") + xlab("COH-1 ChIP signal per 2kb bin")
p3


fountgr<-GRanges(dd[dd$set=="fountain tips +-2kb",])
nonfountgr<-GRanges(dd[dd$set==paste0(">",minDistance,"kb from fountain tips"),])

table(seqnames(fountgr))
table(seqnames(nonfountgr))


## log2 FC ------
salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
salmongr$gene_length<-width(salmongr)

#' Function to assign fountain or non fountain status to gene expression
#'
#' @param salmongr genomic ranges with expression data
#' @param fountgr genomic ranges for fountains
#' @param nonfountgr genomic ranges for non fountain bins
#' @return genomic ranges of TSS with expression data and fountVnonfount status
#' @export
getOverlappingTSS<-function(salmongr,fountgr,nonfountgr){
  salmontss<-resize(salmongr,width=1,fix="start")
  salmontss$fountVnonfount<-NA
  ol<-findOverlaps(salmontss,fountgr)
  salmontss$fountVnonfount[queryHits(ol)]<-"fount"
  fountgr<-fountgr[subjectHits(ol)]
  ol<-findOverlaps(salmontss,nonfountgr)
  salmontss$fountVnonfount[queryHits(ol)]<-"nonfount"
  nonfountgr<-nonfountgr[subjectHits(ol)]
  matched<-list(RNAseq=salmontss,COH1=c(fountgr,nonfountgr))
  return(matched)
}

matched<-getOverlappingTSS(salmongr,fountgr,nonfountgr)
salmontss<-matched[["RNAseq"]]
table(salmontss$fountVnonfount)

COH1ol<-data.frame(matched[["COH1"]])
levels(COH1ol$set)<-c("fountain tips\n +-2kb",paste0(">",minDistance,"kb from\nfountain tips"))
stat.test<-COH1ol  %>% wilcox_test(COH1~set,alternative="greater") %>%
  add_xy_position()
stats_df<-COH1ol %>% group_by(set) %>% summarise(count=n(),COH1=-Inf)

p3b<-ggplot(COH1ol,aes(x=set,y=COH1)) +
  geom_boxplot(outlier.shape=NA,fill="grey90") +
  ylab("COH-1 ChIP signal per 2kb bin") + theme(axis.title.x=element_blank())+
  stat_pvalue_manual(stat.test, label = "p", remove.bracket=F,
                     y=2.3)+
  geom_text(data=stats_df,mapping=aes(label=count),vjust=-0.5)+
  theme(axis.title.x=element_blank())+coord_cartesian(ylim=c(0.5,2.5)) +
  geom_hline(yintercept=median(COH1ol$COH1[COH1ol$set=="fountain tips\n +-2kb"]),
             linetype="dashed")
p3b



toplot<-data.frame(salmontss)[!is.na(data.frame(salmontss)$fountVnonfount),]
toplot$fountVnonfount<-factor(toplot$fountVnonfount)
levels(toplot$fountVnonfount)<-c("fountain tips\n +-2kb",paste0(">",minDistance,"kb from\nfountain tips"))
stats_df<-toplot %>% group_by(fountVnonfount) %>% summarise(count=n(),log2FoldChange=-0.35)
stat.test<-toplot  %>% wilcox_test(log2FoldChange~fountVnonfount, alternative="greater") %>%
  add_xy_position()


p5<-ggplot(toplot,aes(x=fountVnonfount,y=log2FoldChange)) +
  geom_boxplot(fill="grey90",outlier.shape=NA) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_cartesian(ylim=c(-0.35,0.5)) +
  geom_text(data=stats_df,mapping=aes(label=count)) +
  stat_pvalue_manual(stat.test, label = "p",remove.bracket=F,
                     y=0.4,tip.length = 0.01) +
  theme(axis.title.x=element_blank()) +
  ylab(label="log<sub>2</sub>FC")
p5

seed=1243
set.seed(seed)
COH1_wilcox<-list()
LFC_wilcox<-list()
for(i in 1:100){
  # sample matching regions
  dd<-getMatchingSets(fountbins,nonfountbins,metricName="COH1")
  fountgr<-GRanges(dd[dd$set==1,])
  nonfountgr<-GRanges(dd[dd$set==2,])
  # find regions overlapping TSS
  matched<-getOverlappingTSS(salmongr,fountgr,nonfountgr)
  salmontss<-matched[["RNAseq"]]
  COH1ol<-data.frame(matched[["COH1"]])
  # get p-values
  COH1_wilcox[i]<-wilcox.test(COH1ol$COH1[COH1ol$set==1],COH1ol$COH1[COH1ol$set==2],alternative="greater")$p.value
  LFC_wilcox[i]<-wilcox.test(salmontss$log2FoldChange[salmontss$fountVnonfount=="fount"],
                             salmontss$log2FoldChange[salmontss$fountVnonfount=="nonfount"], alternative="greater")$p.value
}


testnames<-c("COH-1 ChIP","RNAseq log<sub>2</sub>FC")
names(testnames)<-c("COH1","LFC")
pvals<-rbind(data.frame(test="COH1",pval=unlist(COH1_wilcox)),
             data.frame(test="LFC",pval=unlist(LFC_wilcox)))
p7<-ggplot(pvals,aes(x=pval)) +geom_histogram()+
  facet_wrap(~test,labeller=labeller(test=testnames)) +
  geom_vline(xintercept=0.05,colour="red")+
  xlab("p-value") +
  theme(strip.text.x = element_markdown())
p7


p<-ggpubr::ggarrange(p2,ggpubr::ggarrange(p3,p3b,ncol=2,widths=c(1.2,1), labels=c("b ","c ")),
                     ggpubr::ggarrange(p5,p7,ncol=2,widths=c(1,1.2),labels=c("d ","e ")),
                     nrow=3,labels=c("a ","",""))
ggsave(paste0("COH1vsRNAseq_matchedCOH1.pdf"),p,device="pdf",width=18,height=21,units="cm")
ggsave(paste0("COH1vsRNAseq_matchedCOH1.png"),p,device="png",width=18,height=21,units="cm")
