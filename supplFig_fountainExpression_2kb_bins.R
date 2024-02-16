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

fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"

rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")



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
ggsave(filename=paste0(finalFigDir,"/supplFig_fountainExpression_2kb_bins.pdf"),
       p, device="pdf",height=29,width=19, units="cm")

#' ######################################
#' ## chromatin states ---------
#' ######################################
#'
#' states<-import.bed(paste0(publicDataDir,"/chromStates_L3_Evans2016_ce11.bed"))
#' #domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
#' #seqlevels(domains)<-seqlevels(states)
#'
#' stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
#'              "#c4d69c","#05ff7f","#ceff65","#fd0082",
#'              "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
#'              "#528dd4","#814006","#b8cde4","#808080",
#'              "#dcdac3","#c4bc97","#938b54","#141414")
#'
#' getStateOLtable<-function(gr,states){
#'   ol<-data.frame(findOverlaps(gr,states,ignore.strand=T,minoverlap=10))
#'   gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
#'   ol$name<-gr$name[ol$queryHits]
#'   ol$XvA<-ifelse(unlist(strsplit(ol$name,":.*$"))=="chrX","X","A")
#'   ol$state<-factor(states$score[ol$subjectHits],levels=1:20)
#'   ol$width<-width(states)[ol$subjectHits]
#'   df<-ol%>%group_by(state,XvA) %>% summarise(stateFrequency=n(),stateWidth=sum(width))
#'   allStates<-data.frame(state=factor(1:20, levels=1:20))
#'   df<-left_join(allStates,df)
#'   df[is.na(df)]<-0
#'   return(df)
#' }
#'
#' #' Binning by distance
#' #'
#' #' Creates bins of a particular size from a vector of distances
#' #' @param distance vector of distances
#' #' @param maxDist The maximum distance bin (greater distances will be pooled together)
#' #' @param binSize Size of bins (in bp, will be converted to kb for labels)
#' #' @param fullTile If using bins completely tiling genome then need to create an
#' #' additional bin for overlapping which is labelled "-", as immediately preceding and
#' #' following bins will be 0 distance too
#' #' @return vector of binned values
#' binByDistance1<-function(distance,maxDist=50000,binSize=2000,fullTile=F){
#'   #levels(binnedDistance)<-seq(binSize,maxDist,by=binSize)/1000
#'   if(fullTile){
#'     binnedDistance<-cut(distance,breaks=seq(0,maxDist-binSize,by=binSize),include.lowest=F)
#'     levels(binnedDistance)<-seq(binSize+binSize,maxDist,by=binSize)/1000
#'     levels(binnedDistance)<-c(seq(binSize+binSize,maxDist,by=binSize)/1000,paste0(">",(maxDist/1000)),0,binSize/1000)
#'     binnedDistance[is.na(binnedDistance)]<-paste0(">",(maxDist/1000))
#'     binnedDistance[distance<binSize/1000]<-binSize/1000
#'     binnedDistance[distance==0]<-0
#'     binnedDistance<-factor(binnedDistance,levels=c(0, binSize/1000, seq(binSize+binSize,maxDist,binSize)/1000, paste0(">",(maxDist/1000))))
#'   } else {
#'     binnedDistance<-cut(distance,breaks=seq(0,maxDist,by=binSize),include.lowest=F)
#'     levels(binnedDistance)<-seq(binSize,maxDist,by=binSize)/1000
#'     levels(binnedDistance)<-c(seq(binSize,maxDist+binSize,by=binSize)/1000,paste0(">",(maxDist/1000)),0)
#'     binnedDistance[is.na(binnedDistance)]<-paste0(">",(maxDist/1000))
#'     binnedDistance[distance==0]<-0
#'     binnedDistance<-relevel(binnedDistance, ref="0")
#'   }
#'   #table(binnedDistance)
#'   return(binnedDistance)
#' }
#'
#' tiles<-tileGenome(seqlengths(Celegans)[1:6],tilewidth=2000,cut.last.tile.in.chrom = T)
#' tiles<-getUpDownstreamFountainData(tiles,fountains)
#'
#' tiles<-makeGRangesFromDataFrame(tiles,keep.extra.columns = T)
#'
#' tiles$binnedDistance<-binByDistance1(tiles$distanceToFountain,maxDist=10000,binSize =2000,fullTile=T)
#' tiles<-flipDownstreamOrientation(tiles)
#' table(tiles$binnedDistance1,tiles$fountainLocation)
#' tiles$AsymmetryTercile<-factor(tiles$AsymmetryTercile,levels=paste0("Asym.Q",1:3))
#'
#' listdf<-list()
#' for(l in levels(tiles$fountainLocation)){
#'   for(a in levels(tiles$AsymmetryTercile)){
#'     for(d in levels(tiles$binnedDistance1)[c(2:6,8:13)]){
#'       df1<-getStateOLtable(gr=tiles[tiles$AsymmetryTercile==a & tiles$binnedDistance1==d],states)
#'       df1$fountainLocation<-l
#'       df1$AsymmetryTercile<-a
#'       df1$binnedDistance1<-d
#'       listdf[[paste0(l,"_",a,"_",d)]]<-df1
#'     }
#'   }
#' }
#'
#'
#' df<-do.call(rbind,listdf)
#' df$binnedDistance1<-factor(df$binnedDistance1,levels=levels(tiles$binnedDistance1)[c(2:6,8:13)])
#'
#' df$XvA<-factor(df$XvA)
#'
#' bin0<-ceiling(length(levels(df$binnedDistance1))/2)
#'
#' boxdf<-df %>% filter(XvA=="A",AsymmetryTercile %in% c("Asym.Q1","Asym.Q3")) %>%
#'   group_by(binnedDistance1,AsymmetryTercile) %>%
#'   mutate(stateWidth=sum(stateWidth),stateFrequency=sum(stateFrequency)) %>%
#'   filter(binnedDistance1=="0") %>% select(-AsymmetryTercile,-state,-fountainLocation) %>% distinct()
#'
#' p3<-ggplot(df[df$XvA=="A" & df$AsymmetryTercile %in% c("Asym.Q1","Asym.Q3"),],
#'            aes(x=binnedDistance1,y=stateWidth/1e6,fill=state)) +
#'   geom_bar(position="stack",stat="identity") +
#'   scale_fill_manual(values=stateClrs) +
#'   xlab("Binned distance (kb)")+
#'   theme(legend.position = "none") +
#'   facet_grid(rows=vars(AsymmetryTercile))+
#'   ggtitle(paste0("Chromatin state width in autosomal fountains")) +ylab("Bin coverage x Mb")
#' p3<-p3 + geom_bar(data=boxdf,aes(x=binnedDistance1, y=stateWidth/1e6,fill=NULL),stat="identity",
#'                   alpha=0,color="black",linewidth=0.5)
#'
#' p3
#'
#'
#'
#'
