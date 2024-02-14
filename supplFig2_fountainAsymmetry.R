library(GenomicRanges)
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

fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
fountains$AsymmetryTercile<-cut(fountains$Asymetry,quantile(fountains$Asymetry,seq(0,1,1/3)),labels=paste0("QAsym",1:3), include.lowest=T)


#rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")

#fileList<-data.frame(sampleName=c("COH1cs"),
#                     filePath=c(rnaSeqFile))

# control bins
#nonFount<-gaps(fountains)
#nonFount<-resize(nonFount,width=2000,fix="center")

#padjVal=0.05

## prepare table
# salmon<-readRDS(file=rnaSeqFile)
# salmon<-salmon[!is.na(salmon$chr),]



#' Assign up and downstream fountains to each gene
#'
#' @param salmongr GRanges object of gene TSSs with expression score (LFC or TPM)
#' @param fountains GRanges object of fountain tips
#' @return data frame of gene expression and closest fountain data (upstream and
#' downstream)
#' @export
getUpDownstreamFountainData<-function(salmongr,fountains){
  originalMcol=ncol(mcols(salmongr))
  ### get fountains downstream of genes
  salmongrNext<-salmongr
  dtnp<-precede(salmongrNext,fountains,ignore.strand=T) # index of fountain that is preceded by gene
  naidx<-is.na(dtnp)
  sum(naidx) # 28
  salmongrNext<-salmongrNext[!naidx]
  salmongrNext$distanceToFountain<-distance(salmongrNext,fountains[dtnp[!naidx]],
                                            ignore.strand=T)
  mcols(salmongrNext)<-cbind(mcols(salmongrNext),mcols(fountains[dtnp[!naidx]]))
  mcols(salmongrNext)$fountainStrand<-strand(fountains[dtnp[!naidx]])

  #fix overlapping fountains that are wrongly assigned
  ol<-findOverlaps(salmongrNext,fountains,ignore.strand=T)
  salmongrNext[queryHits(ol)]$distanceToFountain<- -1
  mcols(salmongrNext)[queryHits(ol),]<-cbind(mcols(salmongrNext)[queryHits(ol),1:(originalMcol+1)],
                                             mcols(fountains)[subjectHits(ol),],
                                             strand(fountains)[subjectHits(ol)])
  #make bedgraph for verification
  # forBG<-salmongrNext
  # mcols(forBG)<-NULL
  # forBG$name<-salmongrNext$fountainName
  # forBG$score<-salmongrNext$score
  # export(forBG,paste0(rnaSeqDir,"/tracks/tpmByDownstreamFountain.bed"))

  # convert to table
  dfNext<-data.frame(salmongrNext)
  dfNext$fountainLocation<-"downstream"

  #### get fountains upstream of genes
  salmongrPrev<-salmongr
  dtnf<-follow(salmongr,fountains,ignore.strand=T) # index of fountain that is followed by gene
  naidx<-is.na(dtnf)
  sum(naidx) #28 #41
  salmongrPrev<-salmongrPrev[!naidx] #28 #41
  salmongrPrev$distanceToFountain<-distance(salmongrPrev,fountains[dtnf[!naidx]],
                                            ignore.strand=T)
  mcols(salmongrPrev)<-cbind(mcols(salmongrPrev),mcols(fountains[dtnf[!naidx]]))
  mcols(salmongrPrev)$fountainStrand<-strand(fountains[dtnf[!naidx]])

  #fix overlapping fountains that are wrongly assigned
  ol<-findOverlaps(salmongrPrev,fountains,ignore.strand=T)
  salmongrPrev[queryHits(ol)]$distanceToFountain<- -1
  mcols(salmongrPrev)[queryHits(ol),]<-cbind(mcols(salmongrPrev)[queryHits(ol),1:(originalMcol+1)],
                                             mcols(fountains)[subjectHits(ol),],
                                             strand(fountains)[subjectHits(ol)])
  #make bedgraph for verification
  # forBG<-salmongrPrev
  # mcols(forBG)<-NULL
  # forBG$name<-salmongrPrev$fountainName
  # forBG$score<-salmongrPrev$score
  # export(forBG,paste0(rnaSeqDir,"/tracks/tpmbyUpstreamFountain.bed"))

  # convert to table
  dfPrev<-data.frame(salmongrPrev)
  dfPrev$fountainLocation<-"upstream"

  #### combine upstream and downstream
  df<-rbind(dfPrev,dfNext)
  df$fountainLocation<-factor(df$fountainLocation,levels=c("upstream","downstream"))
  return(df)
}


#' Binning by distance
#'
#' Creates bins of a particular size from a vector of distances
#' @param distance vector of distances
#' @param maxDist The maximum distance bin (greater distances will be pooled together)
#' @param binSize Size of bins (in bp, will be converted to kb for labels)
#' @param fullTile If using bins completely tiling genome then need to create an
#' additional bin for overlapping which is labelled "-", as immediately preceding and
#' following bins will be 0 distance too
#' @return vector of binned values
binByDistance<-function(distance,maxDist=50000,binSize=2000,fullTile=F){
  binnedDistance<-cut(distance,breaks=seq(0,maxDist,by=binSize),include.lowest=F)
  levels(binnedDistance)<-seq(binSize,maxDist,by=binSize)/1000
  if(fullTile){
    levels(binnedDistance)<-c(seq(binSize,maxDist,by=binSize)/1000,paste0(">",(maxDist/1000)),0,"-")
  } else {
    levels(binnedDistance)<-c(seq(binSize,maxDist,by=binSize)/1000,paste0(">",(maxDist/1000)),0)
  }
  binnedDistance[is.na(binnedDistance)]<-paste0(">",(maxDist/1000))
  binnedDistance[distance==0]<-0
  if(fullTile){
    binnedDistance[distance==-1]<-"-"
  } else {
    binnedDistance[distance==-1]<-0
  }
  binnedDistance<-relevel(binnedDistance, ref="0")
  if(fullTile){
    binnedDistance<-relevel(binnedDistance, ref="-")
  }
  #table(binnedDistance)
  return(binnedDistance)
}

#' Flip orientation of bins for downstream fountains
#'
#' Requires data.frame with "binnedDistance" and "fountainLocation" columns.
#' Will flip the "downstream" fountain level order and add a "-" to label
flipDownstreamOrientation<-function(df){
  originalLevels<-levels(df$binnedDistance)
  df$binnedDistance1<-as.vector(df$binnedDistance)
  df$binnedDistance1[df$fountainLocation=="downstream"]<-paste0("-", df$binnedDistance[df$fountainLocation=="downstream"])
  newLevels<-c(paste0("-",rev(levels(df$binnedDistance))),
               levels(df$binnedDistance))
  df$binnedDistance1<-factor(df$binnedDistance1,
                             levels=newLevels)
  levels(df$binnedDistance1)<-gsub("->",">-",levels(df$binnedDistance1))
  return(df)
}


# salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
# salmongr<-sort(salmongr)
# salmongr$geneLength<-width(salmongr)
# salmongr<-resize(salmongr,width=1,fix="start")
# df<-getUpDownstreamFountainData(salmongr,fountains)
#
# df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=10000,binSize=2000)
# table(df$binnedDistance)
#
# df<-flipDownstreamOrientation(df)
# table(df$binnedDistance1)



### looking at 366 TPM - 10 kb zoom----
tpm366<-import(paste0(rnaSeqDir,"/tracks/PMW366_TPM_avr.bed"))
tpm366$length<-width(tpm366)
tpm366tss<-resize(tpm366,width=1,fix="start")
df<-getUpDownstreamFountainData(tpm366tss,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=10000,binSize=2000)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)


df$score<-NULL
df<-left_join(df,data.frame(tpm366)[,c(6,7)],by=c("name"))

table(df$binnedDistance,df$AsymmetryTercile)
table(df$binnedDistance1)


df<-df[df$score>0,]
#df$AsymmetryTercile<-paste0("QAsym",df$AsymmetryTercile)
ignoreBins<-c(">-10","-0",">10","QAsym2")

selected=!df$binnedDistance1 %in% ignoreBins & !(df$AsymmetryTercile %in% ignoreBins)
countLabels_a<-df[selected,]%>%
  dplyr::group_by(binnedDistance1,AsymmetryTercile,.drop=T) %>%
  dplyr::summarise(count=dplyr::n())

# calculate stats
df$log2TPM<-log2(df$score+1)
stat_df<-df[!(df$AsymmetryTercile %in% ignoreBins),] %>%
  rstatix::group_by(AsymmetryTercile,binnedDistance,drop=T) %>% filter(binnedDistance==2) %>%
  rstatix::t_test(log2TPM~binnedDistance1)  %>% filter(group2 %in% 2:10) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group1 %in% ignoreBins) %>%
  add_xy_position(x="binnedDistance1") %>%
  mutate(xmin=xmin-1,xmax=xmax-2,y.position[order(AsymmetryTercile,y.position)])

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p13<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2TPM),
               fill=rep(c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
                        length(unique(countLabels_a$AsymmetryTercile))),
               outlier.color="grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=14,size=3,colour="blue",
            angle=0) +
  xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+1)")+
  ggtitle("2kb bins of log2 TPM (366) by fountain QAsymetry and location - all genes") +
  facet_grid(AsymmetryTercile~.) + coord_cartesian(ylim=c(0,15))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
                     hide.ns = T, color="purple",bracket.size=0.5)
p13
stat_df$p.adj





p<-ggarrange(p13,p13a,ncol=2)
ggsave(filename=paste0(outPath,"/plots/TPMv",regionType,"Asymmetry_distanceBoxplot_QAsymGrp_10kb_all.pdf"),
       p, device="pdf",width=29,height=19, units="cm")



### looking at coh-1 TPM - 10 kb zoom----
tpmCOH1<-import(paste0(rnaSeqDir,"/tracks/PMW828_TPM_avr.bed"))
#names(mcols(tpmCOH1))<-c("wormbaseID","COH1tpm")

tpmCOH1$length<-width(tpmCOH1)
tpmCOH1tss<-resize(tpmCOH1,width=1,fix="start")
df<-getUpDownstreamFountainData(tpmCOH1tss,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=10000,binSize=2000)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)


df$score<-NULL
df<-left_join(df,data.frame(tpmCOH1)[,c(6,7)],by=c("name"))


table(df$binnedDistance,df$AsymmetryQuintile)

df$AsymmetryQuintile<-paste0("QAsym",df$AsymmetryQuintile)
ignoreBins<-c(">-10","-0",">10","QAsym2","QAsym3","QAsym4")

df<-df[df$score>0,]


selected=!df$binnedDistance1 %in% ignoreBins & !(df$AsymmetryQuintile %in% ignoreBins)
countLabels_a<-df[selected,]%>%
  dplyr::group_by(binnedDistance1,AsymmetryQuintile,.drop=T) %>%
  dplyr::summarise(count=dplyr::n())

# calculate stats
df$log2TPM<-log2(df$score+1)
stat_df<-df[!(df$AsymmetryQuintile %in% ignoreBins),] %>%
  rstatix::group_by(AsymmetryQuintile,binnedDistance,drop=T) %>% filter(binnedDistance==2) %>%
  rstatix::t_test(log2TPM~binnedDistance1)  %>% filter(group2 %in% 2:10) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group1 %in% ignoreBins) %>%
  add_xy_position(x="binnedDistance1") %>%
  mutate(xmin=xmin-1,xmax=xmax-2,y.position[order(AsymmetryQuintile,y.position)])

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p14<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2TPM),
               fill=rep(c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
                        length(unique(countLabels_a$AsymmetryQuintile))),
               outlier.color="grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=14,size=3,colour="blue",
            angle=0) +
  xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+1)")+
  ggtitle("2kb bins of log2 COH-1 TPM by fountain QAsymetry and location - all genes") +
  facet_grid(AsymmetryQuintile~.) + coord_cartesian(ylim=c(0,15))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
                     hide.ns = T, color="purple",bracket.size=0.5)
p14
stat_df


p<-ggarrange(p14,p14a,ncol=2)
ggsave(filename=paste0(outPath,"/plots/COH1TPMv",regionType,"Asymmetry_distanceBoxplot_QAsymGrp_10kb_all.pdf"),
       p, device="pdf",width=29,height=19, units="cm")

