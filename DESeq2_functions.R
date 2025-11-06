###############################################################-
# Create creating tx2gene and metadata GRanges object ---------
###############################################################-

#' Create TxDb from Wormbase annotation gtf
#'
#' Convert wormbase annotation gtf to txdb object for particular
#' genome version. File must be named as in wormbase:
#' paste0("c_elegans.PRJNA13758.",genomeVer,".annotations.gtf")
#' @param genomeDir Directory where gtf is stored.
#' @param genomeVer Version of wormbase to use, e.g. "WS285"
#' @return None. sqlite txdb will be written to genomeDir.
#' @export
makeTxDbsqlite<-function(genomeDir,genomeVer){
    #dir.create(paste0(genomeDir,"/annotations"),recursive=T)
    si<-GenomeInfoDb::seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
    GenomeInfoDb::genome(si)<-genomeVer
    GenomeInfoDb::seqnames(si)<-gsub("M","MtDNA",gsub("chr","",GenomeInfoDb::seqnames(si)))
    wstxdb<-GenomicFeatures::makeTxDbFromGFF(file=paste0(genomeDir,
                                        "/c_elegans.PRJNA13758.",
                                        genomeVer,".annotations.gtf"),
                            format="gtf",organism="Caenorhabditis elegans",
                            chrominfo=si)

    saveDb(wstxdb,paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                         ".annotations.sqlite"))
}

#' Create TxDb from Wormbase annotation gtf but save with ce11 genome
#'
#' Convert wormbase annotation gtf to txdb object for particular
#' genome version. File must be named as in wormbase:
#' paste0("c_elegans.PRJNA13758.",genomeVer,".annotations.gtf")
#' @param genomeDir Directory where gtf is stored.
#' @param genomeVer Version of wormbase to use, e.g. "WS285"
#' @return None. sqlite txdb will be written to genomeDir.
#' @export
makeTxDbsqlite_ce11<-function(genomeDir,genomeVer){
  #dir.create(paste0(genomeDir,"/annotations"),recursive=T)
  gtf<-import(paste0(genomeDir,"/c_elegans.PRJNA13758.",
                     genomeVer,".annotations.gtf"))
  si<-GenomeInfoDb::seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqlevels(gtf)<-seqlevels(si)
  seqinfo(gtf)<-seqinfo(Celegans)
  gtf$transcript_id<-gsub("Transcript:","",gtf$transcript_id)
  gtf$gene_id<-gsub("Gene:","",gtf$gene_id)
  export(gtf,paste0(genomeDir,"/c_elegans.PRJNA13758.",
                genomeVer,"_ce11.annotations.gtf"))
  wstxdb<-GenomicFeatures::makeTxDbFromGFF(file=paste0(genomeDir,
                                                       "/c_elegans.PRJNA13758.",
                                                       genomeVer,"_ce11.annotations.gtf"),
                                           format="gtf",organism="Caenorhabditis elegans",
                                           chrominfo=si)

  saveDb(wstxdb,paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                       "_ce11.annotations.sqlite"))
}



#' Make tx2gene object
#'
#' Use txdb from wormbase annotation to make a tx2gene table required to
#' read in salmon data into DEseq2 for gene level counts.
#' @param genomeDir Directory where gtf is stored.
#' @param genomeVer Version of wormbase to use, e.g. "WS285"
#' @return None. sqlite txdb will be written to genomeDir.
#' @export
getTx2Gene<-function(genomeDir,genomeVer){
  if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer, ".annotations.sqlite"))){
    makeTxDbsqlite(genomeDir, genomeVer)
  }
  # load a txdb of wormbase data and create a tx2gene object
  txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  tx2gene$TXNAME<-gsub("Transcript:","",tx2gene$TXNAME)
  #print("tx2gene columns: "); columns(txdb)
  #print("tx2gene keytypes: "); keytypes(txdb)
  return(tx2gene)
}

#' Create GRanges of gene names
#'
#' Create a GRanges object for protein coding genes from wormbase annotation
#' txdb sqlite file. Requires txdb and geneIDs file from wormbase, both
#' saved into genomeDir. Metadata object uses ucsc format seqnames.
#' @param genomeDir Directory where txdb sqlite and geneIDs files are stored.
#' @param genomeVer Version of wormbase to use, e.g. "WS285"
#' @param outPath where to save the metadata rds object
#' @param format "ce11" or "wb" to choose ucsc or wormbase formatted chromosomes
#' @return Metadata object
#' @export
getMetadataGR<-function(genomeDir,genomeVer,outPath=".",format="ce11"){
  if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer, ".annotations.sqlite"))){
    makeTxDbsqlite(genomeDir, genomeVer)
  }
  # load a txdb of wormbase data and create a tx2gene object
  if(!file.exists(paste0(outPath,"/ce11GeneGR_",genomeVer,".rds"))){
    txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                        ".annotations.sqlite"))
    TxptByGene<-transcriptsBy(txdb, by = "gene")
    print("TxptByGene length: "); length(TxptByGene)

    geneGR<-unlist(range(TxptByGene))
    mcols(geneGR)$wormbase<-names(geneGR)
    genedf<-as.data.frame(geneGR)
    genedf$wormbase<-gsub("^Gene:","",genedf$wormbase)

    geneIDs<-read.delim(gzfile(paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".geneIDs.txt.gz")), header=FALSE, sep=",")
    colnames(geneIDs) <- c("species.ID","WormBase.Gene.ID","Public.Name","Sequence.Name","gene.anotation.status","class")
    # can be ignored: david
    #david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")

    ### combine wormbase gene IDs with gene location ###
    metadata<-inner_join(geneIDs, genedf,by=c("WormBase.Gene.ID"="wormbase")) %>%
      dplyr::select(WormBase.Gene.ID,Public.Name,Sequence.Name,class,seqnames,start, end, strand) %>%
      collect %>% GenomicRanges::GRanges()

    names(mcols(metadata))<-c("wormbaseID","publicID","sequenceID","class")

    table(metadata$class)
    saveRDS(metadata,paste0(outPath,"/wbGeneGR_",genomeVer,".rds"))

    seqlevelsStyle(metadata)<-"ucsc"
    seqinfo(metadata)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
    metadata<-sort(metadata)

    saveRDS(metadata,paste0(outPath,"/ce11GeneGR_",genomeVer,".rds"))
  }
  if(format=="ce11"){
    metadata<-readRDS(paste0(outPath,"/ce11GeneGR_",genomeVer,".rds"))
  } else {
    metadata<-readRDS(paste0(outPath,"/wbGeneGR_",genomeVer,".rds"))
  }
  return(metadata)
}

##########################-
# basic QC -----
##########################-

#' Write dds sample statistics to file
#'
#' Writes basic statistics of the sample counts in the DDS DESeq2 object to
#' a log file paste0(outPath,"/txt/", fileNamePrefix,"all_logfile.txt").
#' @param dds DESeq2 object
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output log.
#' @return Stats written to file
#' @export
writeDDSstatsToFile<-function(dds,outPath=".", fileNamePrefix="salmon_"){
  ## basic sample stats
  sink(file=paste0(outPath,"/txt/", fileNamePrefix,
                        "all_logfile.txt"),
            append=FALSE, type="output")
  statsPerSample<-data.frame(t(apply(counts(dds),2,summary)))
  statsPerSample$totalCounts<-colSums(counts(dds))
  rownames(statsPerSample)<-colData(dds)$sampleName
  colnames(statsPerSample)<-c("min", "Q1", "median", "mean", "Q3", "max","totalCounts")
  statsPerSample$zeros <- apply(counts(dds)==0, 2, sum)
  statsPerSample$percZeros <- round(100*statsPerSample$zeros/nrow(counts(dds)),1)
  print(statsPerSample)
  sink()
}


#' Get long list of colours from base R
#'
#' Get long list of colour from base R so that one can easily plot large numbers
#' of samples. Selecting all colours with "4" in their name, excluding whites and
#' greys. If 'paired' is FALSE, will return list of these colours with the 4
#' stripped off. If paired is TRUE, will return a light and dar version of the
#' colour list.
#' @param paired logical value: should colour list have pairs of colours (default=F)
#' @return a vector of colours
#' @export
getColours<-function(paired=F){
  clrs<-colors()[grep("4$",colors())]
  clrs<-clrs[grep("(grey)|(white)|(gray)",clrs,invert=T)]
  clrs<-gsub("4$","",clrs)
  pairedClrs<-as.vector(sapply(clrs,function(x) {c(paste0(x,"3"),paste0(x,"4"))}))
  if(paired==T){
    return(pairedClrs)
  } else {
    return(clrs)
  }
}


#' Compare raw RNAseq counts across samples
#'
#' Use boxplots to compare raw RNAseq counts across samples, coloured
#' by sample group. DDS object must have column called sampleName with
#' name of samples.
#' @param dds DESeqDataSet object
#' @param groupingVariable Name of variable in colData that indicates the grouping
#' variable of interest to colour samples by group
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Data written to pdf file
#' @export
compareSampleCounts<-function(dds,groupingVariable,outPath=".",fileNamePrefix="salmon_"){
  grpClrs<-getColours()
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                  "QCcounts_Box_Density.pdf"), width=8,height=8,paper="a4")
  ## box plots
  sampleGroups<-colData(dds)[[groupingVariable]]
  epsilon <- 1 # pseudo-count to avoid problems with log(0)
  par(mar=c(5,8,4,2))
  boxplot(log2(counts(dds) + epsilon), col=grpClrs[as.numeric(sampleGroups)],
          pch=".",horizontal=TRUE, cex.axis=1,
          las=1, ylab=NULL, names=colData(dds)$sampleName,
          xlab=paste0("log2(counts ","+",epsilon,")"))
  par(mar=c(5,4,4,2))

  ## density plots
  plotDensity(log2(counts(dds) + epsilon), lty=1, lwd=2, xlab="log2 counts per gene",
              col=grpClrs[as.numeric(sampleGroups)])
  grid()
  legend("topright", legend=colData(dds)$sampleName, lwd=2,cex=0.8,
         col=grpClrs[as.numeric(sampleGroups)])
  dev.off()
}


#' Plot estimations of size factors and dispersion
#'
#' Quality control for DESeq2 model  - check size factors to make sure median
#' ratio normalisation is appropriate (distribution should be "normal"). And
#' check estimated dispersion values are centered around the estimated curve.
#' @param dds DESeqDataSet object
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Data written to pdf file
#' @export
modelQC<-function(dds,outPath=".",fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                  "QCmodel_sizeFactors_dispersion.pdf"), width=8,height=8,paper="a4")
  counts <- counts(dds, normalized = TRUE)
  p<-degCheckFactors(counts) # Distribution of gene ratios used to calculate Size Factors
  print(p)
  # checks median ratio normalization is appropriate - should be "normally" disributed
  plotDispEsts(dds)
  #Dispersion, accounts for a gene’s variance and mean expression level. Calculated by Var = μ + α*μ^2
  #https://hbctraining.github.io/DGE_workshop_salmon/lessons/04_DGE_DESeq2_analysis.html
  dev.off()
}



#' Plot heatmaps of sample to sample and genes by sample clustering
#'
#' Quality control to check how well samples cluster by group. DESeqDataSet object
#' must have a colData column called sampleName with the names of the samples.
#' @param dds DESeqDataSet object
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Data written to pdf file
#' @export
clusteringSamples_heatmap<-function(dds,outPath=".",fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                  "QCsampleClustering_heatmaps.pdf"), width=8,height=8,paper="a4")
  # sample to sample
  vsd <- vst(dds, blind=TRUE)
  idx<-order(rowMeans(counts(dds,normalized=T)),decreasing = T)[1:min(10000,nrow(assay(dds)))]
  colnames(vsd)<-colData(dds)$sampleName
  sampleDists <- dist(t(assay(vsd)[idx,]))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colData(dds)$sampleName
  colnames(sampleDistMatrix) <- colData(dds)$sampleName
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  p<-pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  print(p)
  #samples by top 500 expressed genes
  select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:500]
  df <- data.frame(colData(dds)[,3:ncol(colData(dds))])
  rownames(df)<-colData(vsd)$sampleName
  p<-pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=TRUE, annotation_col=df,
           main="Top 500 expressed genes (avr counts) - clustered by columns")
  print(p)
  dev.off()
}



#' Plot clustering of covariates
#'
#' Quality control to check how well samples cluster by group. DESeqDataSet object
#' must have a colData column called sampleName with the names of the samples.
#' @param dds DESeqDataSet object
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Data written to pdf file
#' @export
checkCovariates<-function(dds,outPath=".",fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                  "QCcovariateClustering.pdf"), width=8,height=8,paper="a4")
  print("plotting qc from DEGreport package")
  counts <- counts(dds, normalized = TRUE)
  design <- as.data.frame(colData(dds))
  colnames(counts)<-design$sampleName
  rownames(design)<-design$sampleName
  degCovariates(log2(counts+0.5),design)
  degCorCov(design)
  dev.off()
}


#' Plot PCA of samples
#'
#' Quality control to check how well samples cluster by group. DESeqDataSet object
#' must have a colData column called sampleName with the names of the samples.
#' @param dds DESeqDataSet object
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Data written to pdf file
#' @export
plotPCAbyCovariates<-function(dds,outPath=".",fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                  "QC_PCA.pdf"), width=8,height=8,paper="a4")
  colourByVar<-colnames(sampleTable)[grep("fileName|sampleName",colnames(sampleTable),invert=T)]
  plotlist<-list()
  rowAvr<-rowMeans(assay(dds))
  idx<-order(rowAvr,decreasing = T)[1:min(10000,nrow(assay(dds)))]
  vsd <- vst(dds[idx], blind=TRUE)
  for(varToUse in colourByVar){
    p1<-plotPCA(vsd, intgroup=c(varToUse),ntop=5000)+ggtitle(varToUse)+
      theme(legend.text = element_text(size=8))
    plotlist[[varToUse]]<-p1
  }
  p<-ggarrange(plotlist=plotlist,ncol=1,nrow=3)
  print(p)
  dev.off()
}


#' Plot sample to sample correlations
#'
#' Quality control to check how well samples cluster by group. DESeqDataSet object
#' must have a colData column called sampleName with the names of the samples.
#' @param dds DESeqDataSet object
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @param groupingVariable Name of variable in colData that indicates the grouping variable
#' of interest for performing sub-group correlations
#' @param plotPDFs boolean: should output be as pdf? Default False produces smaller png files.
#' @return Data written to pdf file
#' @export
plotCorrelations<-function(dds,outPath=".",fileNamePrefix="salmon_",
                           groupingVariable=NULL,plotPDFs=F){
  print("plotting pairwise correlation between genes")
  #Define a function to draw a scatter plot for a pair of variables (samples) with density colors
  plotFun <- function(x,y){
    dns <- densCols(x,y);
    points(x,y, col=dns, pch=".", panel.first=grid());
    #  abline(a=0, b=1, col="brown")
  }

  corFun <- function(x,y){
    par(usr = c(0, 1, 0, 1))
    txt <- as.character(format(cor(x, y), digits=2))
    text(0.5, 0.5, txt, cex = 3*cor(x, y),
         col=colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)[9*(cor(x,y))])
  }
  epsilon<-1
  df<-as.data.frame(log2(counts(dds) + epsilon))
  colnames(df)<-colData(dds)$sampleName
  if(!is.null(groupingVariable)){
    groupsOI<-unique(colData(dds)[[groupingVariable]])
    for (grp in groupsOI){
      if(plotPDFs==T){
        pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.pdf"),
            width=8,height=8,paper="a4")
      } else {
        png(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.png"),
            width=8,height=8,units="in",res=150)
      }
      idx<-colData(dds)[[groupingVariable]] %in% c(grp)
      if(sum(idx)>1){
        pairs(df[,idx], panel=plotFun, lower.panel=corFun, labels=colnames(df)[idx], main=grp)
      }
      dev.off()
    }
  }
  if(plotPDFs==T){
    pdf(file=paste0(outPath,"/plots/",fileNamePrefix, "all_correlations.pdf"),
        width=12,height=12,paper="a4")
  } else {
    png(file=paste0(outPath,"/plots/",fileNamePrefix, "all_correlations.png"),
        width=24,height=24,units="in",res=150)
  }
  #tmp<-data.frame(do.call(rbind,strsplit(colnames(df),"_")))
  #df<-df[,order(tmp$X3,tmp$X2)]
  pairs(df, panel=plotFun, lower.panel=corFun, labels=colnames(df), main="All samples")
  dev.off()
}


#' Get shrunken Log2 Fold change values
#'
#' Shrink Log2 Fold Change values using lfcShrink() function from DESeq2 and add
#' metadata about gene location and names
#' @param dds DESeqDataSet object
#' @param res DESeqResults object
#' @param metadata GRanges object with metadata about alternative gene names
#' @param contrastName String indicating the name of the contrast used to make res
#' @param shrinkMethod Method to use for LFC shrinking
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Data written to csv and rds files and returned as DESeqResults object
#' @export
getShrunkenLFC<-function(dds, res, metadata, contrastName , shrinkMethod="ashr",
                         outPath=".", fileNamePrefix="salmon_"){
  resLFC<-DESeq2::lfcShrink(dds, type=shrinkMethod, res=res)
  #class(resLFC)
  ### add metadata for further processing
  resLFC$wormbaseID<-rownames(resLFC)
  idx<-match(rownames(resLFC),metadata$wormbaseID)
  resLFC$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
  resLFC$start<-as.vector(start(metadata))[idx]
  resLFC$end<-as.vector(end(metadata))[idx]
  resLFC$strand<-as.vector(strand(metadata))[idx]
  resLFC$publicID<-as.vector(metadata$publicID)[idx]
  resLFC$sequenceID<-as.vector(metadata$sequenceID)[idx]
  resLFC$class<-as.vector(metadata$class)[idx]
  if("oscillating" %in% names(mcols(metadata))){
    resLFC$oscillating<-as.vector(metadata$oscillating)[idx]
  }
  #resLFC$entrezID<-as.vector(metadata$entrezID)[idx]
  saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix, contrastName,
                             "_DESeq2_fullResults.rds"))

  sink(file=paste0(outPath,"/txt/",fileNamePrefix, contrastName,
                   "_logfile.txt"), append=TRUE,
       type="output")
  cat(paste0("Number of genes that change expression in ",contrastName," at different padj cutoffs:\n"))
  print(summary(resLFC))
  print(summary(resLFC,alpha=0.05))
  print(summary(resLFC,alpha=0.01))
  sink()

  #resLFC$entrezID<-as.vector(metadata$entrezID)[idx]
  saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix, contrastName,
                             "_DESeq2_fullResults.rds"))

  #export csv with ordered results
  write.csv(resLFC[order(resLFC$padj),],
            file=paste0(outPath,"/txt/", fileNamePrefix,contrastName,
                        "_DESeq2_resultsTable.csv"),
            quote=F,row.names=F)
  return(resLFC)
}

#' Plot QC of pvalues and filtering threshold
#'
#' Plot histogram of pvalues and ajusted pvalues. Expect to see uniform distribution
#' and elevated counts in smallest bin.
#' Also plot criteria used to choose independent filtering threshold.
#' @param resLFC DESeqResults object with shrunken LFC
#' @param contrastName String indicating the name of the contrast used to make res
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plots written to pdf file
#' @export
plotPvalQC<-function(resLFC,contrastName,outPath=".", fileNamePrefix="salmon_"){
  # histograms of pvals
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastName,
                  "_QCdge_pvalHist_filtThresh.pdf"), width=5,height=11,paper="a4")
  par(mfrow=c(3,1))
  if("oscillating" %in% names(resLFC)){
    oscFac<-factor(resLFC$oscillating,levels=c("no","Meeuse","Latorre","Meeuse&Latorre"))
    pvalLevels<-cut(resLFC$pvalue,seq(0,1,0.1),include.lowest=T)
    barplot(table(oscFac, pvalLevels), legend.text=levels(oscFac),
            main=paste0(contrastName," pvalue"))
    pvalLevels<-cut(resLFC$padj,seq(0,1,0.1),include.lowest=T)
    barplot(table(oscFac, pvalLevels), legend.text=levels(oscFac),
            main=paste0(contrastName," adjusted pvalue"))
  } else {
    hist(resLFC$pvalue,breaks=9,main=paste0(contrastName," pvalue"))
    hist(resLFC$padj,breaks=9, main=paste0(contrastName," adjusted pvalue"))
  }
  # filtering threshold
  plot(metadata(resLFC)$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter",main="Threshold for independant filtering of results")
  lines(metadata(resLFC)$lo.fit, col="red")
  abline(v=metadata(resLFC)$filterTheta)
  legend("topright",legend=paste0("Read count \nthreshold: ",
                                  round(metadata(resLFC)$filterThreshold,2)))

  par(mfrow=c(1,1))
  dev.off()
}







#' Plot heatmap of counts of significantly changing genes
#'
#' Plot heatmap of heirarchically clustered samples based on significantly changing
#' genes.
#' @param dds DESeqDataSet object
#' @param resLFC DESeqResults object
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
plotHClust_significant<-function(dds,resLFC,contrastName, padjVal, outPath=".",
                                 fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastName,
                  "_hclust_mostChanged.pdf"), width=8,height=11,paper="a4")
  gene.kept <- rownames(resLFC)[resLFC$padj <= padjVal & !is.na(resLFC$padj) &
                                  abs(resLFC$log2FoldChange)>0]
  epsilon<-1
  if(length(gene.kept)>10){
    # Retrieve the normalized counts for gene of interest
    countTable.kept <- log2(counts(dds, normalized=TRUE) + epsilon)[gene.kept, ]
    dim(countTable.kept)
    colnames(countTable.kept)<-colData(dds)$sampleName

    # Perform the hierarchical clustering with
    # A distance based on Pearson-correlation coefficient
    # and average linkage clustering as agglomeration criteria
    par(cex.main=0.8)
    heatmap.2(as.matrix(countTable.kept),
              scale="row",
              hclust=function(x)hclust(x,method="average"),
              distfun=function(x) as.dist((1-cor(t(x)))/2),
              margin=c(7,0),
              trace="none",
              density="none",
              labRow="",
              #labCol = names(countTable.kept),
              cexCol=1,
              keysize=1,
              main=paste0(contrastName," ",nrow(countTable.kept),
                      " changed genes (p<",padjVal,", |lfc|>0)"))
    par(cex.main=1.2)
  }
  dev.off()
}

#' Plot barplots of counts of 20 most significantly changing genes
#'
#' Plot barplots of 20 most significantly changing genes.
#' @param dds DESeqDataSet object
#' @param resLFC DESeqResults object
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param groupingVariable Name of variable in colData that indicates the grouping variable
#' of interest to colour barplots by group
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
plotCountsMostSig_barplot<-function(dds,resLFC,contrastName,groupingVariable, outPath=".",
                                    fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastName,
                  "_topGenes_normCounts.pdf"), width=8,height=11,paper="a4")
  grpClrs<-getColours()
  par(mfrow=c(5,2))
  selectedGenes <- rownames(resLFC)[order(resLFC$padj)][1:20]
  print("selectedGenes are:")
  print(selectedGenes)
  groupsOI<-colData(dds)[[groupingVariable]]
  if(length(selectedGenes)>1){
    for (g in selectedGenes) {
      g.title <- paste0(g, ":", metadata$publicID[which(metadata$wormbaseID == g)])
      barplot(counts(dds, normalized=TRUE)[g,],
              col=c(grpClrs)[as.numeric(groupsOI)],
              main=g.title, las=2, cex.names=1,
              names.arg=colData(dds)$sampleName)
      # las:2 perpendicular; substring redundant prefix from sample name:8
      legend("topleft",legend=levels(groupsOI),
             fill=c(grpClrs), cex=0.7)
    }
  }
  dev.off()
}


#' Make MA plots
#'
#' Plot expression change vs mean expression (MA plot) for both shurnken and
#' nonshrunken values. A few genes can be highlighted and labelled by passing a
#' list of their WBGene numbers
#' @param res DESeqResults object before shrinking
#' @param resLFC DESeqResults object after shrinking
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @param highlight Vector of WBGene IDs of genes you want to label in MA plot
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
makeMAplots<-function(res,resLFC,contrastName, padjVal, highlight=NULL ,outPath=".",
                      fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastName, "_MAplots_results.pdf"), width=8,height=11,paper="a4")
  par(mfrow=c(2,1))
  plotMA(res, main=paste0(contrastName," uncorrected LFC, threshold=", padjVal), ylim=c(-5,5), alpha=padjVal, colSig = "red3")
  if(!is.null(highlight)){
    highlight_BaseMean <- res[highlight,]$baseMean
    points(highlight_BaseMean, res[highlight,]$log2FoldChange, col="darkblue", pch=1)
    text(highlight_BaseMean, 0.3+res[highlight,]$log2FoldChange, names(highlight), col="black")
  }
  plotMA(resLFC, main=paste0(contrastName," apeglm shrunk LFC, threshold=", padjVal), ylim=c(-5,5), alpha=padjVal, colSig = "red3",cex=1)
  if("oscillating" %in% names(resLFC)){
    points(resLFC[resLFC$oscillating!="no","baseMean"],resLFC[resLFC$oscillating!="no","log2FoldChange"],
           col="grey",pch=16,cex=0.5)
  }
  if(!is.null(highlight)){
    points(highlight_BaseMean, resLFC[highlight,]$log2FoldChange, col="darkblue", pch=1)
    text(highlight_BaseMean, 0.3+resLFC[highlight,]$log2FoldChange, names(highlight), col="black")
  }
  par(mfrow=c(1,1))
  dev.off()
}




#' Plot Log2 Fold Change by chromosome
#'
#' Plot Log2 Fold Change by chromosome
#' @param resLFC DESeqResults object after shrinking
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
plotLFCbyChromosome<-function(resLFC,contrastName, padjVal ,outPath=".",
                              fileNamePrefix="salmon_"){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastName,
                  "_boxPlots_expnByChr.pdf"), width=8,height=5,paper="a4")
  chrName<-factor(resLFC$chr)
  geneCounts<-table(chrName)

  boxplot(log2FoldChange~chrName,data=resLFC,varwidth=TRUE,outline=FALSE,notch=TRUE,
          main=paste0("Expression changes in ", contrastName), ylab="log2 Fold Change",
          col=c(rep("grey",5),"grey"),xlab="chromosome (number of genes)",
          names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
  #stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
  abline(h=0,lty=2,col="blue")
  dev.off()
}


#' Plot counts per chromosome
#'
#' Counts of up and down regulated genes per chromosome
#' @param resLFC DESeqResults object after shrinking
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @param lfcVal Log2 fold change to be used for choosing significant genes
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
plotCountsPerChr<-function(resLFC,contrastName, padjVal, lfcVal,outPath=".",
                           fileNamePrefix="salmon_"){
  summaryByChr<-function(resLFC,padjVal,lfcVal) {
    up<-resLFC[resLFC$padj < padjVal & resLFC$log2FoldChange > lfcVal,]
    down<-resLFC[resLFC$padj < padjVal & resLFC$log2FoldChange < -lfcVal, ]
    allChr<-as.data.frame(rbind(up=table(up$chr),down=table(down$chr)))
    allChr$autosomes<-rowSums(allChr[,1:5])
    allChr$total<-rowSums(allChr[,1:6])
    rownames(allChr)<-paste0(rownames(allChr),"_p",padjVal,"_lfc",lfcVal)
    return(allChr)
  }
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastName,
                  "_barplots_countsUpDownByChr.pdf"), width=8,height=5,paper="a4")
  cnts<-summaryByChr(resLFC,padjVal=padjVal,lfcVal=lfcVal)

  barplot(as.matrix(cnts)[,1:6], beside=T, legend.text=c("up","down"),
  args.legend = list(x = "topleft"),
  main=paste0(contrastName," up/down regulated genes padj<",padjVal," |lfc|>",lfcVal))
  dev.off()
  # get counts for different thresholds
  sink(file=paste0(outPath,"/txt/", fileNamePrefix, contrastName,
                   "_logfile.txt"),append=TRUE, type="output")
  cat("Summary by Chr: \n")
  cat("\np=0.05, LFC=0: \n")
  print(summaryByChr(resLFC,padjVal=0.05,lfcVal=0))
  cat("\np=0.05, LFC=0.5: \n")
  print(summaryByChr(resLFC,padjVal=0.05,lfcVal=0.5))
  cat("\np=0.05, LFC=1: \n")
  print(summaryByChr(resLFC,padjVal=0.05,lfcVal=1))

  cat("\np=0.01, LFC=0: \n")
  print(summaryByChr(resLFC,padjVal=0.01,lfcVal=0))
  cat("\np=0.01, LFC=0.5: \n")
  print(summaryByChr(resLFC,padjVal=0.01,lfcVal=0.5))
  cat("\np=0.01, LFC=1: \n")
  print(summaryByChr(resLFC,padjVal=0.01,lfcVal=1))
  sink()

}




#' Do volcano plot
#'
#' Volcano plot with -log10(pvalue) vs log2(foldChange)
#' @param resLFC DESeqResults object after shrinking
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @param lfcVal Log2 fold change to be used for choosing significant genes
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
makeVolcanoPlot<-function(resLFC,contrastName, padjVal, lfcVal,outPath=".",
                          fileNamePrefix="salmon_"){
  # get point colours
  keyvals<-rep('black', nrow(resLFC))
  names(keyvals)<-rep('NS',nrow(resLFC))
  keyvals[which(resLFC$padj<padjVal & abs(resLFC$log2FoldChange)>lfcVal)]<-'red'
  names(keyvals)[which(resLFC$padj<padjVal & abs(resLFC$log2FoldChange)>lfcVal)]<-paste0('p<',padjVal,' |lfc|>',lfcVal)
  if("oscillating" %in% names(resLFC)){
   # get point shapes
   shapeVals<-rep(19,nrow(resLFC))
   names(shapeVals)<-rep("Not oscillating",nrow(resLFC))
   shapeVals[which(resLFC$oscillating == "Meeuse")]<-3
   names(shapeVals)[which(resLFC$oscillating == "Meeuse")]<-"Meeuse"
   shapeVals[which(resLFC$oscillating == "Latorre")]<-4
   names(shapeVals)[which(resLFC$oscillating == "Latorre")]<-"Latorre"
   shapeVals[which(resLFC$oscillating == "Meeuse&Latorre")]<-8
   names(shapeVals)[which(resLFC$oscillating == "Meeuse&Latorre")]<-"Meeuse&Latorre"
  } else {
    shapeVals<-NULL
  }
  sigUp<-sum(resLFC$padj<padjVal & resLFC$log2FoldChange>lfcVal)
  sigDown<-sum(resLFC$padj<padjVal & resLFC$log2FoldChange< -lfcVal)
  p1<-EnhancedVolcano(resLFC,
                        lab=rownames(resLFC),
                        labSize=0.5,
                        labCol="#11111100",
                        x="log2FoldChange",
                        y="padj",
                        selectLab=rownames(resLFC)[12366],
                        xlim=c(-5.5,5.5),
                        ylim=c(0,65),
                        title= paste0(contrastName),
                        titleLabSize = 12,
                        subtitle=NULL,
                        caption = paste0(nrow(resLFC), ' expressed genes.\n ',sigUp, " up, ",sigDown," down."),
                        captionLabSize = 10,
                        pCutoff=padjVal,
                        FCcutoff=lfcVal,
                        xlab=bquote(~Log[2]~'fold change'~.(contrastName)),
                        ylab=bquote(~-Log[10]~adjusted~italic(P)),
                        #.legend=c('NS','P & Log2 FC'),
                        #legendLabels=c('NS', expression(p-value<padjVal~and~log[2]~FC>1)),
                        legendPosition = 'right',
                        legendLabSize = 8,
                        legendIconSize = 2.0,
                        axisLabSize=10,
                        colCustom=keyvals,
                        #col = c("black", "red"),
                        shapeCustom=shapeVals,
                        colAlpha=0.5,
                        pointSize = 1.0)
    #dev.off()
    if(plotPDFs==T){
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastName,
                             "_volcanoPlot_allGenes.pdf"), plot=p1,
             device="pdf", width=12,height=8,units="cm")
    } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastName,
                             "_volcanoPlot_allGenes.png"), plot=p1,
             device="png", width=12,height=8,units="cm")
    }
}




#' @param resLFC DESeqResults object after shrinking
#' @param contrastName String indicating the name of the contrast used to make resLFC
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @param lfcVal Log2 fold change to be used for choosing significant genes
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plot written to pdf file
#' @export
makeLFCtracks<-function(resLFC,contrastName, padjVal, lfcVal,outPath=".",
                        fileNamePrefix="salmon_"){
  ## make GRanges for LFC
  resGR<-GenomicRanges::GRanges(seqnames=resLFC$chr,
                                IRanges::IRanges(start=resLFC$start,
                                                 end=resLFC$end),
                                strand=resLFC$strand)
  seqlengths(resGR)<-seqlengths(Celegans)[1:6]
  mcols(resGR)<-resLFC[,c("wormbaseID","log2FoldChange","padj")]

  names(mcols(resGR))[names(mcols(resGR))=="log2FoldChange"]<-"score"
  resGR<-sort(resGR,ignore.strand=TRUE)

  #https://github.com/hochwagenlab/hwglabr2/wiki/Apply-function-to-GRanges-scores-in-genomic-tiles
  # make bedgraph
  forBG<-resGR
  mcols(forBG)<-mcols(forBG)[,c("wormbaseID","score")]
  colnames(mcols(forBG))<-c("name","score")
  seqinfo(forBG)<-ce11seqinfo
  export(forBG,paste0(outPath,"/tracks/",fileNamePrefix,contrastName,"_lfc.bedGraph"),
         format="bedGraph")

  # make bigwig
  forBW<-disjoin(forBG,ignore.strand=T)
  oldf<-as.data.frame(findOverlaps(forBW,forBG,ignore.strand=T))
  oldf$scorePerSubBp<-forBG$score[oldf$subjectHits]/width(forBG)[oldf$subjectHits]
  oldf$scorePerQuery<-width(forBW)[oldf$queryHits]*oldf$scorePerSubBp
  score<-oldf %>% group_by(queryHits) %>% summarise(score=mean(scorePerQuery))
  forBW$score<-score$score
  export(forBW,paste0(outPath,"/tracks/",fileNamePrefix,contrastName,"_lfc.bw"),
         format="bigwig")

  ## bed file for significant genes
  idx<-which(resGR$padj<padjVal)
  forBed<-resGR[idx]
  mcols(forBed)<-mcols(forBed)[,c("wormbaseID","score")]
  colnames(mcols(forBed))<-c("name","score")
  seqinfo(forBed)<-ce11seqinfo
  #NaIdx<-is.na(forBed$score)
  #forBed$score[NaIdx]<-0
  export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,contrastName,"_lfc_p",gsub("^0.","",padjVal),".bed"),
         format="bed")

}


#' Remove oscillating genes
#'
#' Remove genes that oscillate during larval development. Based on two studies
#' from Meeuse et al. (2020) and Laorre (2015).
#' @param dds DESeq2 object
#' @param remove One of c("either","both","Meeuse","Latorre"). If "either"
#' removes genes present in either one of these studies. "both" removes genes
#' present in both studies. And "Meeuse" and "Latorre" only remove genes from
#' that one particular study
#' @return DESeq2 object with rows corresponding to specificed genes removed
#' @export
removeOscillating<-function(dds,remove="either"){
  if(!file.exists(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) |
     !file.exists(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))){
    source("./processPublished.R")
  }
  latorre<-read.delim(file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
                      header=T,sep="\t")
  meeuse<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))
  if(remove=="either"){
    idlist=unique(c(latorre$wormbaseID,meeuse$wormbaseID))
  }
  if(remove=="both"){
    idx<-which(meeuse$wormbaseID %in% latorre$wormbaseID)
    idlist<-meeuse$wormbaseID[idx]
  }
  if(remove=="Meeuse"){
    idlist<-meeuse$wormbaseID
  }
  if(remove=="Latorre"){
    idlist<-latorre$wormbaseID
  }
  idx<-rowData(dds)$gene %in% idlist
  dds<-dds[!idx,]
  print(paste0("Removing ",sum(idx)," oscillating genes from ",remove," dataset"))
  return(dds)
}

#' Remove oscillating genes
#'
#' Remove genes that oscillate during larval development. Based on two studies
#' from Meeuse et al. (2020) and Laorre (2015).
#' @param metadata GRanges object
#' @return GRanges object with genes tagged by oscillating info
#' @export
tagOscillating<-function(metadata){
  if(!file.exists(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) |
     !file.exists(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))){
    source("./processPublished.R")
  }
  latorre<-read.delim(file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
                      header=T,sep="\t")
  meeuse<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))
  idlist=unique(c(latorre$wormbaseID,meeuse$wormbaseID))
  df<-data.frame(wormbaseID=unique(c(latorre$wormbaseID,meeuse$wormbaseID)),
                 oscillating=NA)
  idx<-which(df$wormbaseID %in% latorre$wormbaseID)
  df$oscillating[idx]<-"Latorre"
  idx<-which(df$wormbaseID %in% meeuse$wormbaseID)
  df$oscillating[idx]<-"Meeuse"
  idx<-which(df$wormbaseID %in% meeuse$wormbaseID & df$wormbaseID %in% latorre$wormbaseID)
  df$oscillating[idx]<-"Meeuse&Latorre"
  idx1<-which(metadata$wormbaseID %in% df$wormbaseID)
  idx2<-match(metadata$wormbaseID[idx1], df$wormbaseID)
  metadata$oscillating<-"no"
  metadata$oscillating[idx1]<-df$oscillating[idx2]
  return(metadata)
}
