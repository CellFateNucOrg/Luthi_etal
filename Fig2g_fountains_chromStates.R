library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)


projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseeq_DGE")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
#rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")

regionType="fountains"



#######################-
## functions -----------
#######################-

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

  # convert to table
  dfNext<-data.frame(salmongrNext)
  dfNext$fountainLocation<-"downstream"

  #### get fountains upstream of genes
  salmongrPrev<-salmongr
  dtnf<-follow(salmongr,fountains,ignore.strand=T) # index of fountain that is followed by gene
  naidx<-is.na(dtnf)
  sum(naidx) #28
  salmongrPrev<-salmongrPrev[!naidx] #28
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
  #levels(binnedDistance)<-seq(binSize,maxDist,by=binSize)/1000
  if(fullTile){
    binnedDistance<-cut(distance,breaks=seq(0,maxDist-binSize,by=binSize),include.lowest=F)
    levels(binnedDistance)<-seq(binSize+binSize,maxDist,by=binSize)/1000
    levels(binnedDistance)<-c(seq(binSize+binSize,maxDist,by=binSize)/1000,paste0(">",(maxDist/1000)),0,binSize/1000)
    binnedDistance[is.na(binnedDistance)]<-paste0(">",(maxDist/1000))
    binnedDistance[distance<binSize/1000]<-binSize/1000
    binnedDistance[distance==0]<-0
    binnedDistance<-factor(binnedDistance,levels=c(0, binSize/1000, seq(binSize+binSize,maxDist,binSize)/1000, paste0(">",(maxDist/1000))))
  } else {
    binnedDistance<-cut(distance,breaks=seq(0,maxDist,by=binSize),include.lowest=F)
    levels(binnedDistance)<-seq(binSize,maxDist,by=binSize)/1000
    levels(binnedDistance)<-c(seq(binSize,maxDist+binSize,by=binSize)/1000,paste0(">",(maxDist/1000)),0)
    binnedDistance[is.na(binnedDistance)]<-paste0(">",(maxDist/1000))
    binnedDistance[distance==0]<-0
    binnedDistance<-relevel(binnedDistance, ref="0")
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



######################################
########## chromatin states ---------
######################################

states<-import.bed(paste0(projectDir,"/publicData/chromStates_L3_Evans2016_ce11.bed"))
#domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
#seqlevels(domains)<-seqlevels(states)

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")

getStateOLtable<-function(gr,states){
  ol<-data.frame(findOverlaps(gr,states,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$XvA<-ifelse(unlist(strsplit(ol$name,":.*$"))=="chrX","X","A")
  ol$state<-factor(states$score[ol$subjectHits],levels=1:20)
  ol$width<-width(states)[ol$subjectHits]
  df<-ol%>%dplyr::group_by(state,XvA) %>% dplyr::summarise(stateFrequency=dplyr::n(),stateWidth=sum(width))
  allStates<-data.frame(state=factor(1:20, levels=1:20))
  df<-left_join(allStates,df)
  df[is.na(df)]<-0
  return(df)
}

tiles<-tileGenome(seqlengths(Celegans)[1:6],tilewidth=2000,cut.last.tile.in.chrom = T)

tiles<-getUpDownstreamFountainData(tiles,fountains)
tiles<-makeGRangesFromDataFrame(tiles,keep.extra.columns = T)

tiles$binnedDistance<-binByDistance(tiles$distanceToFountain,maxDist=10000,binSize =2000,fullTile=T)
table(tiles$binnedDistance,tiles$fountainLocation)

listdf<-list()
for(l in levels(tiles$fountainLocation)){
  for(a in levels(tiles$AsymmetryQuintile)){
    for(d in levels(tiles$binnedDistance)[1:6]){
      df1<-getStateOLtable(gr=tiles[tiles$fountainLocation==l & tiles$AsymmetryQuintile==a & tiles$binnedDistance==d],states)
      df1$fountainLocation<-l
      df1$AsymmetryQuintile<-a
      df1$binnedDistance<-d
      listdf[[paste0(l,"_",a,"_",d)]]<-df1
    }
  }
}


df<-do.call(rbind,listdf)
df$binnedDistance<-factor(df$binnedDistance,levels=levels(tiles$binnedDistance)[1:6])

df$XvA<-factor(df$XvA)


bin0<-ceiling(length(levels(df$binnedDistance))/2)

boxdf<-df %>% dplyr::filter(XvA=="A") %>%
  dplyr::group_by(binnedDistance) %>%
  dplyr::mutate(stateWidth=sum(stateWidth),stateFrequency=sum(stateFrequency)) %>%
  dplyr::filter(binnedDistance=="0") %>%
  dplyr::select(-AsymmetryQuintile,-state,-fountainLocation) %>%
  dplyr::distinct()

p1<-ggplot(df[df$XvA=="A",],
           aes(x=binnedDistance,y=stateWidth/1e6,fill=state)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_manual(values=stateClrs) +
  theme(legend.position = "none",
        title =element_text(size=5),
        axis.title = element_text(size=10)) +
  ggtitle(paste0("Chromatin state width in autosomal fountains")) +
  ylab(expression(paste("Bins coverage x",10^6," bp"))) +
  xlab(paste0("Fountains (n=",length(fountains),")"))
p1








#### genome wide distribution of chromatin states (autosomes)
tile<-unlist(tileGenome(seqlengths(Celegans),tilewidth = 2000))
ddff<-getStateOLtable(tile,states)
ddff$bin2kb<-factor(0)
p2<-ggplot(ddff[ddff$XvA=="A",],
           aes(x=bin2kb,y=stateWidth/1e7,fill=state)) +
  geom_col(position="stack") +
  scale_fill_manual(values=stateClrs) +
  theme(legend.position = "none",
        title =element_text(size=5),
        axis.title = element_text(size=10)) +
  ggtitle(paste0("ChrA")) +
  ylab(expression(paste("Genome-wide coverage x",10^7," bp"))) +
  xlab("")

p2

p<-ggpubr::ggarrange(p1,p2,ncol=2,widths=c(0.7,0.3))
ggsave(paste0(finalFigDir,"/fig2g_fountains_chromStates.pdf"), p, device="pdf",
       width=8,height=8, units="cm")

