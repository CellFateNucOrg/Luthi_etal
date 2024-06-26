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




getStateOLtable<-function(gr,states){
  ol<-data.frame(findOverlaps(gr,states,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$XvA<-ifelse(unlist(strsplit(ol$name,":.*$"))=="chrX","X","A")
  ol$state<-factor(states$score[ol$subjectHits],levels=1:20)
  ol$width<-width(states)[ol$subjectHits]
  df<-ol%>%group_by(state,XvA) %>% summarise(stateFrequency=n(),stateWidth=sum(width))
  allStates<-data.frame(state=factor(1:20, levels=1:20))
  df<-left_join(allStates,df)
  df[is.na(df)]<-0
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
binByDistance1<-function(distance,maxDist=50000,binSize=2000,fullTile=F){
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


#' Pairwise faceted violin/boxplot with stats
#'
#' Takes a df with columns sampleName, regionType and log2FoldChange and
#' plots violin/boxplot for the two regionTyps faceted by sampleName.
#' t_test is performed on each pair and significance indicated. Total number
#' of regions used also indicated in grey on bottom.
#' @param df with sampleName, regionType and log2FoldChange columns
#' @param gr regions used to plot (needed to calculate width)
#' @param yvar Name of variable to plot on y axis (default is Log2FoldChange)
#' @param ymin bottom of y scale used by coord_cartesian
#' @param ymax top of y scale used by coord_cartesian
#' @param facet_by Name of variable to facet by (default is sampleName). Can
#' be set to NULL to eliminate facetting
#' @param geneSet Name of gene set used in plot title e.g. "all" or "significant"
#' @return plot
#' @export
pairwiseBoxPlotFunc<-function(df,gr,yvar="log2FoldChange",ymin=-1,ymax=1,facet_by="sampleName",geneSet="all"){
  df$sampleName<-factor(df$sampleName, levels=unique(df$sampleName))
  df$regionType<-factor(df$regionType)

  stat_df<-df %>% rstatix::group_by(sampleName,drop=T) %>%
    rstatix::mutate(column=get(yvar)) %>%
    rstatix::wilcox_test(column~regionType) %>%
    rstatix::adjust_pvalue(method="BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x="regionType") %>%
    rstatix::mutate(y.position=0.9*ymax)
  if(!is.null(facet_by)){
    label_df<-df%>% dplyr::group_by(.data[[facet_by]],regionType) %>% dplyr::summarise(count=dplyr::n())
  } else {
    label_df<-df%>% dplyr::group_by(regionType) %>% dplyr::summarise(count=dplyr::n())
  }
  print(stat_df,width=Inf)

  p<-ggplot(df,aes(x=regionType,y=get(yvar))) +
    #geom_violin(width=0.9,mapping=aes(fill=regionType))+
    geom_boxplot(width=0.2,mapping=aes(fill=regionType),outlier.shape=NA)+
    scale_fill_manual(values=c("pink","lightblue"))+
    coord_cartesian(ylim=c(ymin,ymax)) + ylab(yvar) +
    stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,hide.ns = T,
                       color="purple", bracket.size=0.5,bracket.short=0.1,tip.length=0.001)+
    geom_text(data=label_df,aes(label=count),y=ymin,size=3,colour="darkgrey") +
    ggtitle(paste0(yvar," of ",geneSet," genes in fountains (",width(gr[1])/1000,"kb)"))
  if(ymin<0){
    p<-p+geom_hline(yintercept=0,colour="red",linetype=2)
  }
  if(!is.null(facet_by)){
    p<-p+facet_wrap(.~get(facet_by),nrow=3)
  }
  return(p)
}
