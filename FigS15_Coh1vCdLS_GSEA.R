
library("rtracklayer")
library("GenomicRanges")
library("ggplot2")
library(readxl)
library(fgsea)
library(GSEAplot)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)


theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=10),
          axis.title.x=element_text(size=10))
)

projectDir="."
rnaSeqDir=paste0(projectDir,"/RNAseq_DGE")
publicDataDir=paste0(projectDir,"/publicData")
# otherDataDir=paste0(projectDir,"/otherData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

outPath="."
genomeVer="WS285"
genomeDir <- paste0("/Users/semple/Documents/MeisterLab/GenomeVer/",genomeVer)

grp="coh1"
filterOscillating=T

source(paste0(projectDir,"/functions.R"))
source(paste0(projectDir,"/DESeq2_functions.R"))


clrs<-getColours()
plotPDFs <- F

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

ce11seqinfo<-seqinfo(Celegans)

#scriptName="compareGeneLists_HumanToWorm"
#fileNamePrefix=paste0(scriptName,"/")
#makeDirs(outPath,dirNameList=c(paste0(c("plots/"),scriptName)))

padjVal=0.05
lfcVal=0
# DEseq2 results files to use
RNAseqRes<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")

metadata<-getMetadataGR(genomeDir,genomeVer,outPath)
metadata<-tagOscillating(metadata)
#write.table(metadata$wormbaseID[metadata$class=="protein_coding_gene"],
#            paste0(outPath,"/allGeneNames.csv"),row.names=F,col.names=F,quote=F)

salmon<-readRDS(RNAseqRes)



#######################-
## ortholist coverage
#######################-

ortholist<-read_excel(paste0(publicDataDir,"/ortholist_master.xlsx"))
length(unique(ortholist$`WormBase ID`)) #8228
length(unique(ortholist$`HGNC Symbol`)) #11433

names(ortholist)<-make.names(names(ortholist))
numProgs=2
filtOrtholist<- ortholist[ortholist$No..of.Programs>=numProgs,]
dim(filtOrtholist)
filtOrtholist<-filtOrtholist[filtOrtholist$WormBase.ID %in% metadata$wormbaseID[metadata$oscillating =="no"],]
dim(filtOrtholist)
length(unique(filtOrtholist$WormBase.ID))


salmon<-as.data.frame(salmon)[as.data.frame(salmon)$wormbaseID %in% filtOrtholist$WormBase.ID,]


##############-
## Weiss .... Merkenschlager Cornelia de lange paper
##############-
#https://www.nature.com/articles/s41467-021-23141-9
#Neuronal genes deregulated in Cornelia de Lange Syndrome respond to removal and re-expression of cohesin
#Felix D. Weiss, Lesly Calderon, Yi-Fang Wang, Radina Georgieva, Ya Guo, Nevena Cvetesic, Maninder Kaur, Gopuraja Dharmalingam, Ian D. Krantz, Boris Lenhard, Amanda G. Fisher & Matthias Merkenschlager
#Nature Communications volume 12, Article number: 2919 (2021)
options(tibble.width=Inf)
cdl<-read_excel(paste0(publicDataDir,"/Supplementary Data 3. Gene expression from CdLS NeuN+ RNA-seq.xlsx"),skip=3 )

head(cdl)

cdl<-dplyr::inner_join(cdl,filtOrtholist,by=join_by("Gene.Symbol"=="HGNC.Symbol"))

cdl$padj[is.na(cdl$padj)]<-1

cdlUp<-cdl[cdl$padj<0.05 & cdl$log2FoldChange>0,]
cdlDown<-cdl[cdl$padj<0.05 & cdl$log2FoldChange<0,]

sigGenes<-list()
sigGenes[["cdlUp"]]<-cdlUp$WormBase.ID
sigGenes[["cdlDown"]]<-cdlDown$WormBase.ID



lapply(sigGenes,length)

gseaTbl<-list()
leadEdgeTbl<-NULL

# making ranks (ranks go from high to low)
ranks <- salmon$log2FoldChange
names(ranks) <- salmon$wormbaseID
ranks<-sort(ranks,decreasing=T)

fgseaRes <- fgsea(sigGenes, ranks, minSize=5, maxSize = 5000)

fgseaRes[order(padj), ]
#fgseaRes

gseaTbl<-plotGseaTable(sigGenes, ranks, fgseaRes,gseaParam=0.1,render=F,
                              colwidths = c(5, 3, 0.8, 0, 1.2))
pathName<-unlist(fgseaRes[,"pathway"])[2]
pathway<-sigGenes[[pathName]]

stats<-ranks

#' Custom plotEnrichment
#'
#' Based on code for plotEnrichment function from fgsea package
#' @param pathway list of pathways where each item as named by the pathway it
#' contains and has a vector of gene names.
#' @param stats Statistical feature used to sort the genes e.g. LFC, ordered by
#' magnitude
#' @param gseaParam gsea parameter (?)
#' @param ticksize size of tick
#' @return plot
#' @export
plotEnrichment1<-function(pathway, stats, fgseaRes, pathName, gseaParam = 1, ticksSize = 0.2){
  toPlotLFC<-data.frame(x=1:length(stats),y=stats)

  pd<-plotEnrichmentData(pathway, stats, gseaParam)

  g<-with(pd,
       ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="green") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      linewidth=ticksSize,color="darkgreen") +
         geom_hline(yintercept=posES, colour="red", linetype="dashed") +
         geom_hline(yintercept=negES, colour="red", linetype="dashed") +
         geom_hline(yintercept=0, colour="black") +
         theme(
           panel.background = element_blank(),
           panel.grid.major=element_line(color="grey92")
         ) +
         labs(x="rank", y="enrichment score"))

  g<-g+labs(title=paste0("COH-1 cleavage enrichment of ",pathName, " genes")) +
    geom_vline(xintercept=sum(sort(ranks)>0), colour="grey40")+
    annotate("text", x=length(ranks)*0.75, y=fgseaRes[pathway==pathName,ES]*0.9,
             label= paste(paste(paste(c("padj","NES","size"),
                                      fgseaRes[pathway==pathName,c(round(padj,3),round(NES,1),size)],
                                      sep=":"),collapse=", ")))
  g1<-ggplot(toPlotLFC,aes(x=x, y=y)) + geom_bar(stat="identity") +
    xlab("Rank of Log2 ( COH-1cs/TEVcontrol ) High->Low") + ylab("Log2FC")
  p<-ggpubr::ggarrange(g1,g,nrow=2,heights=c(1,2))
  p
}


gseaList<-list()

for (pathName in unlist(fgseaRes[,"pathway"])) {
  gseaList[[pathName]]<-plotEnrichment1(sigGenes[[pathName]], ranks,fgseaRes,pathName) #+
}
pdf(file=paste0(finalFigDir, "/FigS15_Coh1vCdLS_GSEA_",numProgs,"progs.pdf"),
    width=5,height=10,paper="a4")
p<-gridExtra::marrangeGrob(grobs=gseaList, nrow=2, ncol=1,top=paste0("orthologs found with at least ",numProgs," programs"))
print(p)
dev.off()

