library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(InteractionSet)
library(plyranges)
library(GenomicInteractions)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(tidyverse)
library(monocle3)
library(cowplot)
library(readxl)
library(ggpmisc)




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
otherDataDir=paste0(projectDir,"/otherData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

source(paste0(projectDir,"/functions_plotting.R"))



## fountains -----
fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
nonfount<-resize(gaps(fountains),width=2000,fix="center")
export.bed(nonfount,"~/Downloads/nonFount_2kb.bed")
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
seqlevels(fountains)<-seqlevels(Celegans)
fountains<-resize(fountains,width=6000,fix="center")
mcols(fountains)<-NULL
nonfount<-resize(gaps(fountains),width=6000,fix="center")


## RNAseq -------
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")
salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]

# salmon$log2FoldChange_mag<-salmon$log2FoldChange*10
# salmon[!duplicated(salmon$sequenceID),]
# write_delim(data.frame(salmon[,c("sequenceID","log2FoldChange_mag")]),"/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/CelEst/coh1_data.txt", delim="\t")

## Metadata --------
metadata<-readRDS("wbGeneGR_WS285.rds")
seqlevelsStyle(metadata)<-"UCSC"
seqinfo(metadata)<-seqinfo(Celegans)
metadata<-metadata[seqnames(metadata)!="chrM" & !is.na(metadata$wormbaseID)]
metadata$fountDist<-mcols(distanceToNearest(metadata,fountains,ignore.strand=T))$distance


## ATACseq ------
tissueATAC<-import("/Volumes/external.data/MeisterLab/publicData/ATACseq/Durham2021_L2tissueATAC_GSE157017/allTissuesConsensusNamed.bigBed")


fountdf<-data.frame(fountains)
nonfountdf<-data.frame(nonfount)
for(t in unique(tissueATAC$name)){
  fountdf[,paste0(t,"_ATAC")]<-countOverlaps(fountains,tissueATAC[tissueATAC$name==t])
  nonfountdf[,paste0(t,"_ATAC")]<-countOverlaps(nonfount,tissueATAC[tissueATAC$name==t])
}

atacCols<-grep("_ATAC",colnames(fountdf))
sum(is.na(fountdf[,atacCols]))
fountdf$numTissues<-rowSums(fountdf[,atacCols]!=0)
fountdf$totalNumPeaks<-rowSums(fountdf[,atacCols])
fountdf$meanNumPeaksPerTissue<-rowMeans(as.matrix(fountdf[,atacCols]),na.rm=T)
fountdf$sdNumPeaksPerTissue<-rowSds(as.matrix(fountdf[,atacCols]),na.rm=T)
fountdf$cvNumPeaksPerTissue<-fountdf$sdNumPeaksPerTissue/fountdf$meanNumPeaksPerTissue
fountdf$fountVnonfount<-"fountain tip"
dim(fountdf)

atacCols<-grep("_ATAC",colnames(nonfountdf))
sum(is.na(nonfountdf[,atacCols]))
nonfountdf$numTissues<-rowSums(nonfountdf[,atacCols]!=0)
nonfountdf$totalNumPeaks<-rowSums(nonfountdf[,atacCols])
nonfountdf$meanNumPeaksPerTissue<-rowMeans(as.matrix(nonfountdf[,atacCols]),na.rm=T)
nonfountdf$sdNumPeaksPerTissue<-rowSds(as.matrix(nonfountdf[,atacCols]),na.rm=T)
nonfountdf$cvNumPeaksPerTissue<-nonfountdf$sdNumPeaksPerTissue/nonfountdf$meanNumPeaksPerTissue
nonfountdf$fountVnonfount<-"control"
dim(nonfountdf)

toPlot<-rbind(fountdf,nonfountdf)
p1<-ggplot(toPlot,aes(x=numTissues)) +
  geom_histogram()+ facet_wrap(.~fountVnonfount)


cds<-readRDS(paste0(publicDataDir,"/Ghaddar2023_cds_baseline_post_sub.rds"))

if(!file.exists(paste0(publicDataDir,"/Ghaddar2023_markerGenes_n1000.csv"))){
  marker_test_res <- top_markers(cds, group_cells_by="cell_type_group",
                               reference_cells=1000, cores=2,
                               genes_to_test_per_group = 1000)
  write.csv(marker_test_res,paste0(publicDataDir,"/Ghaddar2023_markerGenes_n1000.csv"),row.names=F)
  markers<-left_join(marker_test_res,data.frame(metadata),by=join_by("gene_id"=="wormbaseID"))
  markers<-markers[!is.na(markers$start),]
  markergr<-GRanges(seqnames=markers$seqnames,IRanges(start=markers$start,end=markers$end),strand=markers$strand)
  mcols(markergr)<-markers
  mcols(markergr)[c("seqnames","strand","start","end","width")]<-NULL
  saveRDS(markergr,paste0(publicDataDir,"/Ghaddar2023_markerGenes_n1000.rds"))
}
markergr<-readRDS(paste0(publicDataDir,"/Ghaddar2023_markerGenes_n1000.rds"))
markergr<-markergr[markergr$marker_test_q_value<0.05 & seqnames(markergr)!="chrM" ]

markergr$fountDist<-mcols(distanceToNearest(markergr,fountains,ignore.strand=T))$distance


markerdf<-data.frame(markergr)
markerdf<-markerdf %>% #dplyr::filter(cell_group %in% c("Body wall muscle","Coelomocytes","Germline","Hypodermis","Intestine","Neurons")) %>%
  dplyr::filter(!(cell_group  %in% c("Atypical cells", "Head mesodermal cell", "Unassigned") )) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(count<=1)


tb<-data.frame(table(markerdf$cell_group))
colnames(tb)<-c("tissue","genes")
data.tb <- tibble(x = 100000, y = 0, tbl = list(tb))

ggplot(markerdf, aes(x=fountDist,color=cell_group)) +
  stat_ecdf() +
  geom_table(data=data.tb,aes(x=x,y=y,label=tbl)) +
  coord_cartesian(xlim=c(0,150000))

# pseudo bulk
colData(cds)$cell_type_group
rd<-rowData(cds)
rd<-left_join(data.frame(rd),data.frame(metadata),by=join_by("id"=="wormbaseID"))
rd<-rd[!is.na(rd$start),]
rdgr<-GRanges(seqnames=rd$seqnames,IRanges(start=rd$start,end=rd$end),strand=rd$strand)
mcols(rdgr)<-rd
mcols(rdgr)[c("seqnames","strand","start","end","width")]<-NULL
rdgr<-rdgr[seqnames(rdgr)!="chrM" ]
rdgr$fountDist<-mcols(distanceToNearest(rdgr,fountains,ignore.strand=T))$distance

ggplot(data.frame(rdgr),aes(x=fountDist,num_cells_expressed)) + geom_bin2d(bins=100)
marker_res1 <- top_markers(cds, group_cells_by="cell_type_group",
                           reference_cells=10000, cores=2, marker_sig_test=FALSE,
                           genes_to_test_per_group = 10000)

marker_res1 %>% group_by(gene_id) %>% mutate(numTissues=n())




# Cao 2017 data -----------------
caoData<-"/Users/semple/Documents/MeisterLab/Datasets/Cao_2017_Science_scRNAseq_L2/aam8940_cao_sm_tables_s1_to_s14.xlsx"


maxTissue<-read_excel(caoData,sheet=6,skip=1)
maxCellType<-read_excel(caoData,sheet=7,skip=1)
maxNeuronCluster<-read_excel(caoData,sheet=8,skip=1)

## maxTissue
maxTissue<-inner_join(maxTissue,data.frame(metadata),by=join_by("gene_id"=="wormbaseID"))
mt<-maxTissue[as.numeric(maxTissue$qval)<0.05,]

tb<-mt %>% group_by(max.tissue) %>% summarise(genes=n())
data.tb <- tibble(x = 100000, y = 0, tbl = list(tb))
ggplot(mt, aes(x=fountDist,color=max.tissue)) +
  stat_ecdf() +
  geom_table(data=data.tb,aes(x=x,y=y,label=tbl)) +
  coord_cartesian(xlim=c(0,150000))

## maxCellType
maxCellType<-inner_join(maxCellType,data.frame(metadata),by=join_by("gene_id"=="wormbaseID"))
mt<-maxCellType[as.numeric(maxCellType$qval)<0.05,]

tb<-mt %>% group_by(max.cell.type) %>% summarise(genes=n())
data.tb <- tibble(x = 100000, y = 0, tbl = list(tb))
ggplot(mt, aes(x=max.cell.type, y=fountDist)) +
  geom_boxplot(notch=T) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=median(mt$fountDist)) +
  coord_cartesian(ylim=c(0,150000)) +
  geom_text(mapping=aes(x=max.cell.type,y=-200,label=genes),data=tb)


## maxNeuronCluster
maxNeuronCluster<-inner_join(maxNeuronCluster,data.frame(metadata),by=join_by("gene_id"=="wormbaseID"))
mt<-maxNeuronCluster[as.numeric(maxNeuronCluster$qval)<0.05,]

tb<-mt %>% group_by(max.cluster) %>% summarise(genes=n())
data.tb <- tibble(x = 100000, y = 0, tbl = list(tb))
ggplot(mt, aes(x=max.cluster, y=fountDist)) +
  ggbeeswarm::geom_beeswarm() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=median(mt$fountDist)) +
  coord_cartesian(ylim=c(0,150000)) +
  geom_text(mapping=aes(x=max.cluster,y=-200,label=genes),data=tb)


## number of tissues
tissueConsensusExpr<-read_excel(caoData,sheet=3,skip=1)
#cellTypeConsensusExpr<-read_excel(caoData,sheet=4)
#neuronClusterConsensusExpr<-read_excel(caoData,sheet=5)

exprThreshold=1
tissues<-colnames(tissueConsensusExpr)[3:ncol(tissueConsensusExpr)]
tissueConsensusExpr$numTissues<-rowSums(tissueConsensusExpr[,tissues]>exprThreshold)
tissueConsensusExpr<-inner_join(tissueConsensusExpr,data.frame(metadata),by=join_by("gene_id"=="wormbaseID"))
tissueConsensusExpr$expr_sd<-apply(tissueConsensusExpr[,tissues],1,sd)

tissueConsensusExpr$numTissues<-factor(tissueConsensusExpr$numTissues)

tb<-tissueConsensusExpr %>% group_by(numTissues) %>% summarise(genes=n())
data.tb <- tibble(x = 100000, y = 0, tbl = list(tb))
ggplot(tissueConsensusExpr, aes(x=fountDist,color=numTissues)) +
  stat_ecdf() +
  geom_table(data=data.tb,aes(x=x,y=y,label=tbl)) +
  coord_cartesian(xlim=c(0,150000)) +
  ggtitle(paste0("Expression threshold = ",exprThreshold))




