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
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")

fileList<-data.frame(sampleName=c("COH1cs"),
                     filePath=c(rnaSeqFile))

# control bins
nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=2000,fix="center")

padjVal=0.05




#########################-
# Daugherty ATAC enhancers -----
#########################-

#' Function for adding count data to plots
count_data <- function (y,ymax=5){
  df <- data.frame(y = ymax, label = length(y))
  return(df)
}


#Genome Res. 2017 Dec;27(12):2096-2107.
# doi: 10.1101/gr.226233.117. Epub 2017 Nov 15.
# Chromatin accessibility dynamics reveal novel functional enhancers in C. elegans
# Aaron C Daugherty 1, Robin W Yeo 1, Jason D Buenrostro 1, William J Greenleaf 1 2, Anshul Kundaje 1 3, Anne Brunet 1 4
# https://pubmed.ncbi.nlm.nih.gov/29141961/
# Suplementary table S3 with added nearest gene and distance to gene

# boxplot of LFC by enhancer type
daughertyEnh<-readRDS(paste0(publicDataDir,"/daugherty2017_L3enhancers_ce11.rds"))


### make gene-centered table -----
daughertyEnhByGene<-daughertyEnh %>% as_tibble() %>%
  dplyr::select(L3_chromHMMState,nearestGene,distanceToNearest) %>%
  dplyr::group_by(nearestGene) %>%
  dplyr::summarise(enhancerTypes=paste0(unique(L3_chromHMMState),collapse=","),
                   enhancerCount=dplyr::n(),
                   enhancerUniqCount=length(unique(L3_chromHMMState)),
                   closestEnhancer=min(distanceToNearest),
                   furthestEnhancer=max(distanceToNearest))
dim(daughertyEnhByGene)


daughertyEnhByGene$type<-daughertyEnhByGene$enhancerTypes
daughertyEnhByGene[daughertyEnhByGene$enhancerUniqCount!="1" &
                     grepl("L3_activeEnhancer", daughertyEnhByGene$enhancerTypes), "type"]<-"L3_mixedActive&Repressed"
daughertyEnhByGene[daughertyEnhByGene$enhancerUniqCount!="1" &
                     !grepl("L3_activeEnhancer", daughertyEnhByGene$enhancerTypes),"type"]<-"L3_mixedRepressed"
table(daughertyEnhByGene$type)
colnames(daughertyEnhByGene)[colnames(daughertyEnhByGene)=="nearestGene"]<-"wormbaseID"
options(pillar.width = Inf)
head(daughertyEnhByGene)



salmon<-readRDS(file=rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
salmongr<-sort(salmongr)
tbl<-left_join(as_tibble(salmongr),daughertyEnhByGene,by="wormbaseID")
tbl$type[is.na(tbl$type)]<-"NoEnhancer"
tmpgr<-resize(width=1,fix="start",makeGRangesFromDataFrame(tbl[tbl$type=="NoEnhancer",]))
dtn<-distanceToNearest(tmpgr,daughertyEnh)
tbl[tbl$type=="NoEnhancer","closestEnhancer"]<-mcols(dtn)$distance
tbl$XvA<-ifelse(tbl$seqnames=="chrX","chrX","autosomes")
tbl$upVdown<-ifelse(tbl$log2FoldChange > 0,"up","down")

tbl$sampleName<-fileList$sampleName
tbl$type<-factor(tbl$type,levels=c("NoEnhancer",
                                           paste0("L3_",c("activeEnhancer",
                                                          "mixedActive&Repressed",
                                                          "repressedEnhancer",
                                                          "H3K27me3Repressed",
                                                          "mixedRepressed"))))
table(tbl$type)

boxlabels<-data.frame(type=factor(levels(tbl$type),levels=c("NoEnhancer",
                                                            paste0("L3_",c("activeEnhancer",
                                                                           "mixedActive&Repressed",
                                                                           "repressedEnhancer",
                                                                           "H3K27me3Repressed",
                                                                           "mixedRepressed")))),
                                                            log2FoldChange=-0.4,
                      labels=c("No enhancer", "Active \nenhancer", "Mixed active \n& repressed",
                               "Repressed \nenhancer","H3K27me3 \nenhancer","Mixed \nrepressed"))


#all genes
p1<-ggplot(tbl,aes(x=type,y=log2FoldChange,fill=type)) +
  geom_boxplot(width=0.5,outlier.colour=NA,alpha=0.7,notch=T) +
  #ggtitle("LFC of genes closest to Daugherty et al. (2017) L3 enhancers") +
  theme(axis.text.x = element_blank(),#element_text(angle=90, vjust=0, hjust=3),
        legend.position = "none") +
  stat_summary(fun.data = count_data,fun.args=c(ymax=0.4), geom = "text",
               position = position_dodge(1), size=3, colour="grey30") +
  geom_hline(yintercept=0,col="grey10",linetype="dashed") +
  coord_cartesian(ylim=c(-0.4,0.4)) +
  scale_fill_manual(values=c("darkgrey","red","purple","blue","darkgreen","#146C94")) +
  ylab(label="Log<sub>2</sub>FC") +
  xlab(label="Closest enhancer types to TSS") +
  geom_text(data=boxlabels,aes(x=type,y=log2FoldChange,label=labels),angle=90,
            hjust=0,vjust=0,size=3)
p1


##################-
# distance to fountains of sig genes ------
##################-
grp="COH1cs"
salmon<-readRDS(file=rnaSeqFile)
salmon<-GRanges(salmon[!is.na(salmon$chr),])
d<-distanceToNearest(salmon,fountains,ignore.strand=T)
salmon$distanceToFount<-mcols(d)$distance
#subsetByOverlaps(salmon,fountains)

salmon$upVdown<-"ns"
salmon$upVdown[salmon$padj<padjVal & salmon$log2FoldChange>0]<-"Up"
salmon$upVdown[salmon$padj<padjVal & salmon$log2FoldChange<0]<-"Down"

salmon<-data.frame(salmon)


stat_box_data <- function(y,upper_limit=0) {
  return(
    data.frame(
      y = upper_limit,
      label = paste(length(y))
    )
  )
}



# Statistical test
stat.test <- salmon %>%
  wilcox_test(distanceToFount ~ upVdown) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "upVdown")
stat.test$p.format <- p_format(
  stat.test$p.adj, accuracy = 1e-7,
  leading.zero = T
)

p2<-ggplot(salmon,aes(x=upVdown,y=distanceToFount)) +
  geom_boxplot(varwidth = F,notch=T,fill="lightgrey",outlier.shape = NA) +
  geom_hline(yintercept=median(salmon$distanceToFount[salmon$upVdown=="ns"]), colour="red",linetype=2)+
  ylab(label="Distance to fountain") + xlab("") +
  coord_cartesian(ylim=c(0,90000))+
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.9,
    size = 3
  )

p2<-p2 + stat_pvalue_manual(stat.test, label = "p.adj",y.position=c(70000,75000,80000),
                            tip.length=0.005,color="purple",size=3)

p2




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

# all genes 2kb
stl<-resultsByGRoverlap(fileList,fountains)
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"Fountain\ntip"

stl<-resultsByGRoverlap(fileList,nonFount)
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"Control"

df<-rbind(fountRes,nonFountRes)

p3<-pairwiseBoxPlotFunc(df,fountains,ymin=-0.8,ymax=0.8,facet_by=NULL)
p3<-p3+xlab("") +ylab(label="Log<sub>2</sub>FC")




#######################-
# get gene specific data -----
#######################-
# gtf downloaded from wormbase
gtf<-rtracklayer::import("/Users/semple/Documents/MeisterLab/GenomeVer/WS285/c_elegans.PRJNA13758.WS285.annotations.gtf")

numTranscripts<-as.data.frame(gtf) %>% dplyr::group_by(gene_id) %>%
  dplyr::summarise(numTranscripts=dplyr::n_distinct(transcript_id))

lengthTranscripts<-as.data.frame(gtf) %>% dplyr::group_by(gene_id) %>%
  filter(type=="transcript") %>%
  dplyr::summarise(maxTxLength=max(width)/1e3,minTxLength=min(width)/1e3,avrTxLength=mean(width)/1e3)

numStartEnd<-as.data.frame(gtf) %>% dplyr::group_by(gene_id) %>%
  filter(type=="transcript")%>%
  dplyr::summarise(numStarts=dplyr::n_distinct(start),
                   numEnds=dplyr::n_distinct(end))

numExons<-gtf %>%  plyranges::filter(type=="exon") %>% unique() %>% plyranges::group_by(gene_id) %>%
  plyranges::summarise(numExons=plyranges::n())

geneData<-dplyr::left_join(numTranscripts,lengthTranscripts,by="gene_id")
geneData<-dplyr::left_join(geneData,numStartEnd,by="gene_id")
geneData<-dplyr::left_join(geneData,as.data.frame(numExons),by="gene_id")
geneData$gene_id<-gsub("Gene:","",geneData$gene_id)
names(geneData)[names(geneData)=="gene_id"]<-"wormbaseID"
geneData

rl<-getListOfResults(fileList,padjVal=0.05)
lapply(rl,dim)
rlg<-lapply(rl,dplyr::left_join,geneData,by="wormbaseID")
lapply(rlg,dim)
res<-do.call(rbind,rlg)



##########################-
# gene features in and out of fountains -----
##########################-

stl<-resultsByGRoverlap(fileList,fountains)
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
fountRes<-do.call(rbind,stl)
fountRes$regionType<-"Fountain\ntip"

stl<-resultsByGRoverlap(fileList,nonFount)
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
nonFountRes<-do.call(rbind,stl)
nonFountRes$regionType<-"Control"

df<-rbind(fountRes,nonFountRes)


gFeatures<-c("numTranscripts","numStarts","numEnds","numExons","avrTxLength","maxTxLength")
pretty_gFeatures<-data.frame(gFeat=gFeatures,
                             pretty=c("Number of transcripts per gene","Number of TSS", "Number of TES",
                    "Number of exons", "Average transcript length (kb)", "Maximum transcript length (kb)"))
plotList<-list()
for(gFeat in gFeatures){
  res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
  ymax=quantile(res[,gFeat],0.999)
  ymin=0
  plotList[[gFeat]]<-pairwiseBoxPlotFunc(df[df$sampleName=="COH1cs",],
                                         fountains,yvar=gFeat, ymin=ymin,
                                         ymax=ymax,facet_by=NULL) +
    ylab(pretty_gFeatures$pretty[pretty_gFeatures$gFeat==gFeat]) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1))
}


#ml<-marrangeGrob(plotList,nrow=1,ncol=6)





# skn-1 DTE ------

res<-data.frame(readRDS(paste0(rnaTxSeqDir,"/sleuth_coh1cs_DTE.RDS")))
skn1<-res$target_id[res$publicID=="skn-1"]
tmp<-res[res$publicID=="skn-1",]
tmp<-tmp[complete.cases(tmp),]
tmp$txt<-round(tmp$qval,2)
tmp$txt[tmp$txt>0.05]<-"ns"
tmp$txtPos<-0.15
tmp$label<-c("a","b")

p5<-ggplot(tmp,aes(x=label,y=b)) +
  geom_bar(stat="identity",fill="lightgrey",color="black") +
  xlab("*skn-1* isoform") + ylab(label="Log<sub>2</sub>FC") +
  geom_errorbar(aes(ymin=b-se_b,ymax=b+se_b),width=0.2) +
  geom_hline(yintercept=0) + coord_cartesian(ylim=c(-0.23,0.16)) +
  geom_text(aes(y=txtPos,label=txt),size=3,color="purple")
p5


# assemble final figure -----------

pnull<-ggplot()

p<-ggpubr::ggarrange(ggpubr::ggarrange(p1,p2,p3,ncol=3,widths=c(0.5,0.25,0.25), labels=c("a ", "b ", "d ")),
                     ggpubr::ggarrange(plotlist=plotList, labels=c("e "),nrow=1,ncol=6),
                     ggpubr::ggarrange(pnull,p5,nrow=1,ncol=2,widths=c(0.9,0.2),labels=c(" ","h ")),nrow=3,
                     heights=c(1,0.9,0.8))
p<-annotate_figure(p, top = text_grob("Isiaka et al., Figure 4", size = 14))
ggsave(paste0(finalFigDir,"/fig4_RNAseq.pdf"), p, device="pdf",
       width=19,height=19, units="cm")


