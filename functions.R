#' Make directories
#'
#' @param path String with path to where the directories should be made
#' @param dirNameList Vector of strings with names of directories to create (can include multilevel directories)
#' @return Creates the directories listed in dirNameList
#' @examples
#' makeDirs(path=".",dirNameList=c("/txt","/rds/sample1"))
#' @export
makeDirs<-function(path,dirNameList=c()) {
  sub("\\/$","",path) #remove directory slash if present in path string
  for (d in dirNameList) {
    if (!dir.exists(paste0(path,"/",d))){  # for alignments
      dir.create(paste0(path,"/",d), recursive=TRUE, showWarnings=FALSE)
    }
  }
}


#' Filter DESeq2 table results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param direction Whether to find genes that are less than (lt), or greater than (gt) the log fold change threshold, or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome ("chrX"), or only autosomes ("autosomes")
#' @param outPath Path to working directory
#' @param filenamePrefix Text to add to filename
#' @return filtered table of results which is also automatically written to disk
#' @export
filterResults<-function(resultsTable, padj=0.05, lfc=0, direction="both",
                        chr="all", outPath=".", filenamePrefix="",
                        writeTable=T) {
  sigGenes<-getSignificantGenes(resultsTable,padj,lfc,direction=direction,chr=chr)
  idx<-resultsTable$wormbaseID %in% sigGenes$wormbaseID
  filtTable<-resultsTable[idx,c("baseMean","log2FoldChange","padj",
                                "wormbaseID","chr","start","end","strand")]
  if(writeTable){
    if(!dir.exists(paste0(outPath,"/txt"))){
      dir.create(paste0(outPath,"/txt"))
    }
    write.csv(filtTable,file=paste0(outPath,"/txt/filtResults_p",
                                    padj,"_",direction,"-","lfc",
                                    lfc,"_",chr,".csv"), row.names=F,
              quote=F)
  }
  return(filtTable)
}


#' Get significant genes from  RNAseq results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param namePadjCol Name of column with adjusted P values
#' @param nameLFCcol Name of column with log fold change values
#' @param direction Whether to find genes that are less than (lt), or greater than (gt) the log fold change threshold, or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome ("chrX"), or only autosomes ("autosomes")
#' @param nameChrCol Name of column with chromosome names.
#' @return Filtered table of significant genes at a certain log fold change and adjusted p value.
#' @export
getSignificantGenes<-function(resultsTable, padj=0.05, lfc=0, namePadjCol="padj",
                              nameLfcCol="log2FoldChange", direction="both",
                              chr="all", nameChrCol="chr", outPath="."){
  if(direction=="both") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & abs(resultsTable[,nameLfcCol])>lfc
  } else if(direction=="gt") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & resultsTable[,nameLfcCol]>lfc
  } else if(direction=="lt") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & resultsTable[,nameLfcCol]<lfc
  } else {
    print("direction must be 'both' to get both tails, \n'gt' to get lfc larger than a specific value, \nor 'lt' to get lfc less than a certain value")
  }
  if(chr=="all"){
    idx<-idx
  } else if(chr=="chrX"){
    idx<-idx & !is.na(resultsTable[,nameChrCol]) & resultsTable[,nameChrCol]=="chrX"
  } else if(chr=="autosomes"){
    idx<-idx & !is.na(resultsTable[,nameChrCol]) & resultsTable[,nameChrCol]!="chrX"
  } else {
    print("chr must be one of 'all', 'chrX' or 'autosomes'")
  }
  filtTable<-resultsTable[idx,]
  return(filtTable)
}



#' Get list of results tables
#'
#' Given a table with columns filePath and sampleName, where filePath contains
#' the full path to a DESeq2 results path, read in all the results files in to a
#' list. if padjVal is not null, the files will also be filtered by significance.
#' Additional filtering can be performed on the Log2FoldChange, direction of
#' Log2FoldChange and X vs autosomes. A sampleName column is added to the results
#' tables to identify which dataset they come from.
#' @param fileList data.frame with columns filePath and sampleName
#' @param padjVal padjusted value used for filtering. if NULL (default) then no filtering.
#' @param lfcVal Log2 fold change value used for filtering (default is 0)
#' @param direction direction of log2 fold change. One of "both" (default), "gt" (greater than),
#' or "lt" (less than).
#' @param chr Chromosomes by which to filter results. Default is "all", otherwise "chrX" or "autosomes".
#' @return List of data.frames with DESeq2 results
#' @export
getListOfResults<-function(fileList,padjVal=NULL,lfcVal=0,direction="both",chr="all"){
  sigTables<-list()
  for (i in 1:nrow(fileList)){
    salmon<-readRDS(fileList$filePath[i])
    salmon$sampleName=fileList$sampleName[i]
    if(!is.null(padjVal)){
      sigTables[[fileList$sampleName[i]]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                                             namePadjCol="padj",
                                                                             nameLfcCol="log2FoldChange",
                                                                             direction="both",
                                                                             chr="all", nameChrCol="chr"))
    } else {
      sigTables[[fileList$sampleName[i]]]<-as.data.frame(salmon)
    }
  }
  return(sigTables)
}


#' Get list of results tables subset by overlaps to a GRanges object
#'
#' Given a table with columns filePath and sampleName, where filePath contains
#' the full path to a DESeq2 results path, read in all the results files, subset by
#' their overlap to the GRanges object. if padjVal is not null, the files will also be filtered by significance.
#' Additional filtering can be performed on the Log2FoldChange, direction of
#' Log2FoldChange and X vs autosomes. A sampleName column is added to the results
#' tables to identify which dataset they come from.
#' @param fileList data.frame with columns filePath and sampleName
#' @param gr GenomicRanges object used to select genes which overlap it
#' @param padjVal Adjusted p value used for filtering. if NULL (default) then no filtering.
#' @param lfcVal Log2 fold change value used for filtering (default is 0)
#' @param direction direction of log2 fold change. One of "both" (default), "gt" (greater than),
#' or "lt" (less than).
#' @param chr Chromosomes by which to filter results. Default is "all", otherwise "chrX" or "autosomes".
#' @param resizeTo Defualt is NULL, no resizing of gr. Alternatively can use "start" ,"center", or "end"
#' to make 1bp gr fixed at those positoins.
#' @return List of GRanges with DESeq2 results
#' @export
resultsByGRoverlap<-function(fileList,gr,padjVal=NULL,lfcVal=0,direction="both",chr="all",
                             resizeTo=NULL){
  sigTables<-list()
  #i=1
  for (i in 1:nrow(fileList)){
    salmon<-readRDS(fileList$filePath[i])
    salmon<-salmon[!is.na(salmon$chr),]
    salmon$sampleName=fileList$sampleName[i]
    if(!is.null(padjVal)){
      salmon<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                         namePadjCol="padj",
                                         nameLfcCol="log2FoldChange",
                                         direction="both",
                                         chr="all", nameChrCol="chr")
    }
    resGR<-GRanges(seqnames=salmon$chr,ranges=IRanges(start=salmon$start,end=salmon$end),
                   strand=salmon$strand)
    idx<-colnames(salmon) %in% c("start","end","strand")
    mcols(resGR)<-salmon[,!idx]
    if(!is.null(resizeTo)){
      resGR<-resize(resGR,width=1,fix=resizeTo)
    }
    sigTables[[fileList$sampleName[i]]]<-subsetByOverlaps(resGR,gr,type="any",ignore.strand=T)
  }
  return(sigTables)
}


#' Pairwise faceted violin/boxplot with stats
#'
#' Takes a df with columns sampleName, regionType and log2FoldChange and
#' plots violin/boxplot for the two regionTyps faceted by sampleName.
#' wilcoxon_test is performed on each pair and significance indicated. Total number
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
pairwiseBoxPlotFunc<-function(df,gr,yvar="log2FoldChange",ymin=-1,ymax=1,facet_by="sampleName",geneSet="all",
                              withViolin=T){
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
  if(withViolin){
    p<-ggplot(df,aes(x=regionType,y=get(yvar))) +
      geom_violin(width=0.9,mapping=aes(fill=regionType))+
      geom_boxplot(width=0.2,mapping=aes(fill=regionType),outlier.shape=NA)
  } else {
    p<-ggplot(df,aes(x=regionType,y=get(yvar))) +
      geom_boxplot(width=0.8,mapping=aes(fill=regionType),outlier.shape=NA)
  }
  p<-p+scale_fill_manual(values=c("pink","lightblue"))+
    coord_cartesian(ylim=c(ymin,ymax)) + ylab(yvar) +
    stat_pvalue_manual(stat_df, label = "p.adj", remove.bracket=F,hide.ns = T,
                       color="purple", bracket.size=0.5,bracket.short=0.1,tip.length=0.001,
                       size=3)+
    geom_text(data=label_df,aes(label=count),y=ymin,size=3,colour="grey30") +
    #ggtitle(paste0(yvar," of ",geneSet," genes in fountains (",width(gr[1])/1000,"kb)"))+
    theme(legend.position="none")
  if(ymin<0){
    p<-p+geom_hline(yintercept=0,colour="red",linetype=2)
  }
  if(!is.null(facet_by)){
    p<-p+facet_wrap(.~get(facet_by),nrow=3)
  }
  return(p)
}
