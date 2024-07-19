library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(InteractionSet)
library(plyranges)
library(GenomicInteractions)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(tidyverse)


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

# for large files from HiC fragment mapping
bigDataDir="/Volumes/external.data/MeisterLab/mdas/enhancer_fragment_mapping"

## Parameters ----
maxDistance=30000
maxKNN=10
fragLibs="828v366"
# data from supl table S3 in paper.
totalTEVlib<-sum(258355543, 163424616)
totalCOH1lib<-sum(208298545, 130155406)
#totalSCC1lib<-sum(317973695, 250444362)
#totalSCC1COH1lib<-sum(172186520, 148507656)
#enhancerSet="jaenes"
enhancerSet="daugherty"



## Enhancers -----

if(enhancerSet=="daugherty"){
  enhancers<-readRDS(paste0(publicDataDir,"/daugherty2017_L3Enhancers_ce11.rds"))
  enhancers<-enhancers[enhancers$L3_chromHMMState=="L3_activeEnhancer",]
}
if(enhancerSet=="jaenes"){
  enhancers<-import(paste0(publicDataDir,"/Jaenes2018_L3activeEnhancers_ce11.bed"))
  enhancers$name<-paste0("peak",1:length(enhancers))
}



## Transcripts ----
gtfurl="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS19.canonical_geneset.gtf.gz"
gtffile<-gsub("\\.gz","",basename(gtfurl))

if(!file.exists(paste0(bigDataDir,"/",gtffile))){
  download.file(gtfurl,dest=basename(gtfurl))
  system(paste0("gunzip ",basename(gtfurl)))
}

if(!file.exists(paste0(bigDataDir,"/transcriptPromoters_100upstream100downstream.RDS"))){
  gtf<-import(gtffile)
  seqlevelsStyle(gtf)<-"UCSC"
  seqlevels(gtf)
  table(gtf$type)
  genes<-sort(gtf[gtf$type=="gene" & gtf$gene_biotype=="protein_coding"])
  txpt<-sort(gtf[gtf$type=="transcript" & gtf$gene_biotype=="protein_coding"])
  toBed<-txpt
  toBed$name<-txpt$transcript_id
  toBed$score<-1
  export(toBed,"transcripts_WBPS19.bed")
  txptPromoters<-promoters(txpt,upstream=100,downstream=100)
  saveRDS(txptPromoters,"transcriptPromoters_100upstream100downstream.RDS")
}
txptPromoters<-readRDS(paste0(bigDataDir,"/transcriptPromoters_100upstream100downstream.RDS"))



## HiC fragments ----
fragments<-import(paste0(bigDataDir,"/HicFrag_ce11.bed"))
strand(fragments)<-"*"



## fountains -----
fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
seqlevels(fountains)<-seqlevels(Celegans)
fountains<-resize(fountains,width=6000,fix="center")



## RNAseq -------
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")
salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
#sleuth<-readRDS("/Users/semple/Documents/MeisterLab/papers/Isiaka_etal/RNAseq_DTE/sleuth_coh1cs_DTE.RDS")



## functions ##########################

#' Function to remove trans and very distant interactions
#' and make a GInteractions object with counts
makeCisGI<-function(pairCount,fragments,maxDist=maxDistance){
  first<-sapply(strsplit(pairCount$V1,"\\|"),"[[",1)
  second<-sapply(strsplit(pairCount$V1,"\\|"),"[[",2)


  idx1<-match(first,fragments$name)
  idx2<-match(second,fragments$name)

  gi<-GInteractions(fragments[idx1],fragments[idx2])
  gi$count<-pairCount$V2
  cisgi<-gi[seqnames(anchors(gi,"first")) == seqnames(anchors(gi,"second"))]
  cisgi<-swapAnchors(cisgi)
  cisgi<-sort(cisgi)
  cisgi<-cisgi[pairdist(cisgi)<maxDist]
  return(cisgi)
}


#' Function to combine counts from two GInteraction objects
combineGIcounts<-function(gi1,gi2){
  ol<-findOverlaps(gi1,gi2)
  # sum counts on shared ranges
  gi1$count[queryHits(ol)]<-gi1$count[queryHits(ol)]+gi2$count[subjectHits(ol)]
  # combined shared an unshared regions
  giall<-c(gi1[queryHits(ol)], gi1[! (gi1 %over% gi2) ], gi2[! (gi2 %over% gi1) ])
  giall<-sort(giall)
  return(giall)
}


#' function to find nearest k neighbours irrespective of orientation
#' and return GInteractions ofbject
findKNN<-function(gr1,gr2,maxDistance=maxDistance,maxKNN=4){
  gr1w<-resize(gr1,width=maxDistance,fix="center")
  ol<-data.frame(findOverlaps(gr1w,gr2,ignore.strand=T))
  ol$distance<-GenomicRanges::distance(gr1[ol$queryHits],gr2[ol$subjectHits],ignore.strand=T)
  ol<-ol %>% dplyr::filter(distance!=0) %>%
    dplyr::group_by(queryHits) %>%
    dplyr::mutate(distanceRank=floor(dplyr::row_number(distance))) %>%
    dplyr::top_n(-maxKNN,wt=distanceRank) %>% dplyr::arrange(queryHits,distanceRank)
  ol
  gi<-GInteractions(gr1[ol$queryHits],gr2[ol$subjectHits])
  gi$distance<-ol$distance
  gi$distanceRank<-ol$distanceRank
  return(gi)
}

###########################



## get control counts ----
if(!file.exists(paste0(bigDataDir,"/GI_366_fragment_pair_counts_",maxDistance/1000,"kb.rds"))){
  pairCount1<-read.delim(paste0(bigDataDir,"/366_3_fragment_pair_counts_",maxDistance/1000,"kb.txt"),sep=" ",header=F)
  counts1<-sum(pairCount1$V2)
  counts1
  pairCount2<-read.delim(paste0(bigDataDir,"/366_4_fragment_pair_counts_",maxDistance/1000,"kb.txt"),sep=" ",header=F)  #163424616
  counts2<-sum(pairCount2$V2)
  counts2
  gi1<-makeCisGI(pairCount1,fragments,maxDistance)
  gi2<-makeCisGI(pairCount2,fragments,maxDistance)
  ctrlgi<-combineGIcounts(gi1,gi2)
  #forBedpe<-ctrlgi
  #forBedpe$name<-paste0(ctrlgi$anchor1.name,"|",ctrlgi$anchor2.name)
  #forBedpe$score<-ctrlgi$count
  #export(forBedpe,paste0("GI_366_fragment_pair_counts_",maxDistance/1000,"kb.bedpe"))
  saveRDS(ctrlgi,paste0(bigDataDir,"/GI_366_fragment_pair_counts_",maxDistance/1000,"kb.rds"))
}


## get coh-1 counts ----
if(!file.exists(paste0(bigDataDir,"/GI_828_fragment_pair_counts_",maxDistance/1000,"kb.rds"))){
  pairCount1<-read.delim(paste0(bigDataDir,"/828_1_fragment_pair_counts_",maxDistance/1000,"kb.txt"),sep=" ",header=F)
  counts1<-sum(pairCount1$V2)
  counts1
  pairCount2<-read.delim(paste0(bigDatDir,"/828_2_fragment_pair_counts_",maxDistance/1000,"kb.txt"),sep=" ",header=F)
  counts2<-sum(pairCount2$V2)
  counts2
  gi1<-makeCisGI(pairCount1,fragments,maxDistance)
  gi2<-makeCisGI(pairCount2,fragments,maxDistance)
  coh1gi<-combineGIcounts(gi1,gi2)
  #forBedpe<-coh1gi
  #forBedpe$name<-paste0(coh1gi$anchor1.name,"|",coh1gi$anchor2.name)
  #forBedpe$score<-coh1gi$count
  #export(forBedpe,paste0("GI_828_fragment_pair_counts_",maxDistance/1000,"kb.bedpe"))
  saveRDS(coh1gi,paste0(bigDataDir,"/GI_828_fragment_pair_counts_",maxDistance/1000,"kb.rds"))
}



## join promoter fragments ------
if(!file.exists(paste0(otherDataDir,"/joinedPromoterFragments_100up100down.bed"))){
  promFrag<-join_overlap_inner(fragments,txptPromoters)
  mcols(promFrag)<-mcols(promFrag)[,c("name","gene_id","transcript_id")]
  # add significance category
  promFrag$upVdown<-"NS"
  upGenes<-salmon$wormbaseID[salmon$padj<0.05 & salmon$log2FoldChange>0]
  downGenes<-salmon$wormbaseID[salmon$padj<0.05 & salmon$log2FoldChange<0]
  promFrag$upVdown[promFrag$gene_id %in% upGenes]<-"up"
  promFrag$upVdown[promFrag$gene_id %in% downGenes]<-"down"
  idx<-match(promFrag$gene_id,salmon$wormbaseID)
  promFrag$LFC<-NA
  promFrag$padj<-NA
  promFrag$LFC[!is.na(idx)]<-salmon$log2FoldChange[na.omit(idx)]
  promFrag$padj[!is.na(idx)]<-salmon$padj[na.omit(idx)]

  # add DTE info
  #idx<-match(promFrag$transcript_id,sleuth$target_id)
  #promFrag$DTE_b[!is.na(idx)]<-sleuth$b[na.omit(idx)]
  #promFrag$DTE_qval[!is.na(idx)]<-sleuth$qval[na.omit(idx)]

  # join fragments associated with the same promoter
  promFrag<-data.frame(promFrag) %>% dplyr::group_by(seqnames,transcript_id) %>%
    dplyr::summarise(start=min(start),end=max(end),strand="*",
              gene_id=paste0(unique(gene_id),collapse=";"),
              fragNames=paste0(unique(name),collapse=";"),
              upVdown=paste0(unique(upVdown),collapse=";"),
              LFC=mean(LFC),
              padj=mean(padj,na.rm=T))%>% as_granges()
#              DTE_b=mean(DTE_b,na.rm=T),
#              DTE_qval=mean(DTE_qval)) %>% as_granges()
  table(promFrag$upVdown) # no mixed categories

  # eliminate duplicate fragment names by joining the relevant promoter names (same fragment
  # is promoter for two different genes)
  promFrag<-data.frame(promFrag) %>% dplyr::group_by(seqnames,fragNames) %>%
    summarise(start=min(start),end=max(end),strand="*",
              gene_id=paste0(unique(gene_id),collapse=";"),
              transcript_id=paste0(unique(transcript_id),collapse=";"),
              upVdown=paste0(unique(upVdown),collapse=";"),
              LFC=mean(LFC), padj=mean(padj)) %>% as_granges()
              #DTE_b=mean(DTE_b,na.rm=T),
              #DTE_qval=mean(DTE_qval,na.rm=T)) %>% as_granges()
  table(promFrag$upVdown)
  # get rid of merged up/down and NS fragments
  promFrag$upVdown<-gsub(";NS|NS;","",promFrag$upVdown)
  #salmon[salmon$wormbaseID %in% c("WBGene00269391","WBGene00006975"),]
  promFrag$upVdown<-gsub("down;up","up",promFrag$upVdown) # only one and upreglated gene was more significant
  table(promFrag$upVdown)
  promFrag<-sort(promFrag)
  forBed<-promFrag
  forBed$name<-promFrag$fragNames
  export(forBed,paste0(otherDataDir,"/joinedPromoterFragments_100up100down.bed"))
  saveRDS(promFrag,file=paste0(otherDataDir,"/joinedPromoterFragments_100up100down.RDS"))
}
promFrag<-readRDS(paste0(otherDataDir,"/joinedPromoterFragments_100up100down.RDS"))


## join enhancer fragments ----
if(!file.exists(paste0(otherDataDir,"/joined",enhancerSet,"EnhancerFragments.bed"))){
  enhFrag<-join_overlap_inner(fragments,enhancers)
  colnames(mcols(enhFrag))<-c("fragName","fragScore","enhName","enhScore")
  # join fragments associated with the same enhancer
  enhFrag<-data.frame(enhFrag) %>% dplyr::group_by(seqnames,enhName) %>%
    dplyr::summarise(start=min(start),end=max(end),strand="*",
              fragNames=paste0(unique(fragName),collapse=";"))  %>% as_granges()

  # eliminate duplicate fragment names by joining the relevant enhancer names (same
  # fragment belongs to more than one enhancer)
  enhFrag<-data.frame(enhFrag) %>% dplyr::group_by(seqnames,fragNames) %>%
    dplyr::summarise(start=min(start),end=max(end),strand="*",
              enhName=paste0(unique(enhName),collapse=";")) %>% as_granges()

  enhFrag<-sort(enhFrag)
  forBed<-enhFrag
  forBed$name<-enhFrag$fragNames
  export(forBed,paste0(otherDataDir,"/joined",enhancerSet,"EnhancerFragments.bed"))
  saveRDS(enhFrag,paste0(otherDataDir,"/joined",enhancerSet,"EnhancerFragments.RDS"))
}
enhFrag<-readRDS(paste0(otherDataDir,"/joined",enhancerSet,"EnhancerFragments.RDS"))


##  find enhancer promoter pairs -----

fgi<-findKNN(enhFrag,promFrag,maxDistance,maxKNN)
table(fgi$distanceRank)
fgi$P_E_distanceMid<-pairdist(fgi)



## add count data -----

if(!exists("ctrlgi")){
  ctrlgi<-readRDS(paste0(bigDataDir,"/GI_366_fragment_pair_counts_",maxDistance/1000,"kb.rds"))
  ctrltotalFrag<-sum(ctrlgi$count)
}

if(!exists("coh1gi")){
  coh1gi<-readRDS(paste0(bigDataDir,"/GI_828_fragment_pair_counts_",maxDistance/1000,"kb.rds"))
  coh1totalFrag<-sum(coh1gi$count)
}

ol<-data.frame(findOverlaps(fgi,ctrlgi))
ol$ctrlCounts<-ctrlgi$count[ol$subjectHits]
ol<-ol %>% dplyr::group_by(queryHits) %>% dplyr::mutate(sumCounts=sum(ctrlCounts))
fgi$ctrlCounts<-0
fgi$ctrlCounts[ol$queryHits]<-ol$sumCounts

ol<-data.frame(findOverlaps(fgi,coh1gi))
ol$coh1Counts<-coh1gi$count[ol$subjectHits]
ol<-ol %>% dplyr::group_by(queryHits) %>% dplyr::mutate(sumCounts=sum(coh1Counts))
fgi$coh1Counts<-0
fgi$coh1Counts[ol$queryHits]<-ol$sumCounts


fgi$fountVnonfount<-"nonfount"
ol<-findOverlaps(anchors(fgi)$first,resize(fountains,width=60000,fix="center"))
fgi$fountVnonfount[queryHits(ol)]<-"fountain"
#ol<-findOverlaps(anchors(fgi)$first,fountains)
#fgi$fountVnonfount[queryHits(ol)]<-"tip"
if(!file.exists(paste0(otherDataDir,"/",enhancerSet,"Enh_Prom_",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,".RDS"))){
  saveRDS(fgi,file=paste0(otherDataDir,"/",enhancerSet,"Enh_Prom_",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,".RDS"))
  export.bedpe(fgi, fn = paste0(otherDataDir,"/",enhancerSet,"Enh_Prom_366_",maxDistance/1000,"kb_maxK",maxKNN,".bedpe"), score = "ctrlCounts")
  export.bedpe(fgi, fn = paste0(otherDataDir,"/",enhancerSet,"Enh_Prom_828_",maxDistance/1000,"kb_maxK",maxKNN,".bedpe"), score = "coh1Counts")
}


## convert to table and analyse ------

fgi<-readRDS(paste0(otherDataDir,"/",enhancerSet,"Enh_Prom_",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,".RDS"))
pec<-data.frame(fgi) #daugherty 22941
pec$contact<-paste0(pec$anchor1.fragNames,"|",pec$anchor2.fragNames)
sum(pec$ctrlCounts==0)
sum(pec$coh1Counts==0)
### get ctrl ranks and remove self contacts
pec<-pec %>% dplyr::group_by(anchor1.fragNames) %>%
  dplyr::mutate(ctrlRank=dplyr::min_rank(dplyr::desc(ctrlCounts)),
                coh1Rank=dplyr::min_rank(dplyr::desc(coh1Counts))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(pec$anchor1.fragNames!=pec$anchor2.fragNames)
table(pec$ctrlRank,pec$coh1Rank)

#normalise counts to library size
pec$ctrlCPM<-pec$ctrlCounts*1e6/totalTEVlib
pec$coh1CPM<-pec$coh1Counts*1e6/totalCOH1lib


pec<-pec %>%  dplyr::filter(ctrlCounts>=10, coh1Counts>=10) %>%
  dplyr::mutate(ratioCoh1Ctrl=coh1CPM/ctrlCPM,
         logRatio=log2((coh1CPM)/(ctrlCPM)))
dim(pec) # daugherty 12806
print(pec,width=Inf)



### contact rank -----
pec$ctrlRank<-factor(pec$ctrlRank,levels=1:maxKNN)
stats_df<-pec %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(),ratioCoh1Ctrl=0)
stats_df
stat.test<-pec  %>% wilcox_test(ratioCoh1Ctrl~ctrlRank,ref.group="1",
                                alternative="two.sided",
                                p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p1<-ggplot(pec,aes(x=ctrlRank,y=ratioCoh1Ctrl)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  geom_signif(y_position=rep(c(2,2.1), nrow(stat.test))[1:nrow(stat.test)],
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0)+
  coord_cartesian(ylim=c(0,2.2)) +
  geom_hline(yintercept=1,color="red") +
  ylab("COH-1 / TEV counts ratio") +
  xlab("Promoters ranked by contact probability in TEV control")
p1



### distance rank-----
pec$distanceRank<-factor(pec$distanceRank,levels=1:maxKNN)
stats_df<-pec %>% dplyr::group_by(distanceRank) %>%
  dplyr::summarise(count=dplyr::n(),ratioCoh1Ctrl=0)
stats_df
stat.test<-pec  %>% wilcox_test(ratioCoh1Ctrl~distanceRank,ref.group="1",
                                alternative="two.sided",
                                p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p2<-ggplot(pec,aes(x=distanceRank,y=ratioCoh1Ctrl)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_hline(yintercept=1,color="red")+
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  geom_signif(y_position=rep(c(1.8,1.9), nrow(stat.test))[1:nrow(stat.test)],
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0)+
  coord_cartesian(ylim=c(0,2))  +
  ylab("COH-1 / TEV counts ratio") +
  xlab("Promoters ranked by distance from enhancer")
p2

### LFC by rank -----
stats_df<-pec %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(),anchor2.LFC=-0.4)
stats_df
stat.test<-pec %>% wilcox_test(anchor2.LFC~ctrlRank,ref.group="1",
                               alternative="two.sided",
                               p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p2a<-ggplot(pec,aes(x=ctrlRank,y=anchor2.LFC)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_hline(yintercept=0,color="red")+
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  geom_signif(y_position=rep(c(0.37,0.39), nrow(stat.test))[1:nrow(stat.test)],
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0)+
  coord_cartesian(ylim=c(-0.4,0.4))  +
  ylab("RNAseq Log<sub>2</sub>FC") +
  xlab("Promoters ranked by contact probability in TEV control")
p2a


df1 = pec[!duplicated(pec$ctrlRank) ,]
p3<-ggplot(pec,aes(x=distance/1000,group=ctrlRank)) + geom_histogram() +
  facet_wrap(.~ctrlRank,ncol=5) +
  xlab("Distance of promoter from enhancer (kb)")+
  ylab("Number of promoters ") +
  ggtitle("Contact probability rank in TEV control") +
  geom_text(data=df1,x=7.5,y=+Inf,aes(label=ctrlRank),
            vjust = "inward", hjust = "inward") +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
p3


stats_df<-pec %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(), ctrlCounts=0)
stats_df
p4<-ggplot(pec,aes(y=ctrlCounts,x=ctrlRank)) + geom_boxplot(outlier.shape=NA) +
  xlab("Promoter contact probability rank in TEV control")+
  ylab("Number of contacts ") +
  coord_cartesian(ylim=c(0,ifelse(enhancerSet=="daugherty",200,150))) +
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
                                            size=2.5)
p4

p<-ggpubr::ggarrange(p1,p2,p2a,p3,p4,nrow=5,ncol=1,heights=c(1,1,1,1.1,0.9),labels=c("a","b","c","d","e"))
ggsave(paste0(finalFigDir,"/",enhancerSet,"Enhancer_Promoter_contactRatio",fragLibs,"_maxK",maxKNN,".pdf"),
       p,device="pdf",height=30, width=13, units="cm")



## stratify by fountain -------
fountLabels=c("In fountain (60kb)","Not in fountain")
names(fountLabels)<-c("fountain","nonfount")
### contact rank -----
pec$ctrlRank<-factor(pec$ctrlRank,levels=1:maxKNN)
stats_df<-pec %>% dplyr::group_by(fountVnonfount,ctrlRank) %>% dplyr::summarise(count=dplyr::n(),ratioCoh1Ctrl=0)
stats_df
stat.test<-pec  %>% dplyr::group_by(fountVnonfount) %>%
  wilcox_test(ratioCoh1Ctrl~ctrlRank,ref.group="1",alternative="two.sided",p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test$y_position=rep(c(1.9,2.1), nrow(stat.test))[1:nrow(stat.test)]

p1<-ggplot(pec,aes(x=ctrlRank,y=ratioCoh1Ctrl)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  ggpubr::geom_signif(data=stat.test,
              aes(xmin=group2,xmax=group2, annotations=p.pretty,
                  y_position=y_position),
              parse=T, size=0, textsize=2.5, tip_length=0, manual=T) +
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  coord_cartesian(ylim=c(0,2.2)) +
  geom_hline(yintercept=1,color="red") +
  ylab("COH-1 / TEV counts ratio") +
  xlab("Promoters ranked by contact probability in TEV control")+
  facet_wrap(.~fountVnonfount, nrow=2, labeller=as_labeller(fountLabels),
             strip.position="right")
p1



### distance rank-----
pec$distanceRank<-factor(pec$distanceRank,levels=1:maxKNN)
stats_df<-pec %>% dplyr::group_by(fountVnonfount, distanceRank) %>%
  dplyr::summarise(count=dplyr::n(),ratioCoh1Ctrl=0)
stats_df
stat.test<-pec  %>% dplyr::group_by(fountVnonfount) %>%
  wilcox_test(ratioCoh1Ctrl~distanceRank,ref.group="1",alternative="two.sided", p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")

stat.test$y_position=rep(c(1.8,1.9), nrow(stat.test))[1:nrow(stat.test)]

p2<-ggplot(pec,aes(x=distanceRank,y=ratioCoh1Ctrl)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_hline(yintercept=1,color="red")+
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  ggpubr::geom_signif(data=stat.test,
                        aes(xmin=group2,xmax=group2, annotations=p.pretty,
                            y_position=y_position),
                        parse=T, size=0, textsize=2.5, tip_length=0, manual=T) +
  coord_cartesian(ylim=c(0,2))  +
  ylab("COH-1 / TEV counts ratio") +
  xlab("Promoters ranked by distance from enhancer") +
  facet_wrap(.~fountVnonfount, nrow=2, labeller=as_labeller(fountLabels),
             strip.position="right")
p2



### LFC by rank -----
stats_df<-pec %>% dplyr::group_by(fountVnonfount, ctrlRank) %>% dplyr::summarise(count=dplyr::n(),anchor2.LFC=-0.4)
stats_df
stat.test<-pec %>% dplyr::group_by(fountVnonfount) %>%
  wilcox_test(anchor2.LFC~ctrlRank,ref.group="1",alternative="two.sided", p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test$y_position=rep(c(0.35,0.38), nrow(stat.test))[1:nrow(stat.test)]

p2a<-ggplot(pec,aes(x=ctrlRank,y=anchor2.LFC)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_hline(yintercept=0,color="red")+
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  ggpubr::geom_signif(data=stat.test,
                        aes(xmin=group2,xmax=group2, annotations=p.pretty,
                            y_position=y_position),
                        parse=T, size=0, textsize=2.5, tip_length=0, manual=T) +
  coord_cartesian(ylim=c(-0.4,0.4))  +
  ylab("RNAseq Log<sub>2</sub>FC") +
  xlab("Promoters ranked by contact probability in TEV control")+
  facet_wrap(.~fountVnonfount, nrow=2, labeller=as_labeller(fountLabels),
             strip.position="right")
p2a




df1 = pec[!duplicated(pec$ctrlRank) ,]
p3<-ggplot(pec,aes(x=distance/1000,group=ctrlRank)) + geom_histogram() +
  facet_wrap(.~ctrlRank,ncol=5) +
  xlab("Distance of promoter from enhancer (kb)")+
  ylab("Number of promoters ") +
  ggtitle("Contact probability rank in TEV control") +
  geom_text(data=df1,x=7.5,y=+Inf,aes(label=ctrlRank),
          vjust = "inward", hjust = "inward") +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
p3


stats_df<-pec %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(), ctrlCounts=0)
stats_df
p4<-ggplot(pec,aes(y=ctrlCounts,x=ctrlRank)) + geom_boxplot(outlier.shape=NA) +
  xlab("Promoter contact probability rank in TEV control")+
  ylab("Number of contacts ") +
  coord_cartesian(ylim=c(0,ifelse(enhancerSet=="daugherty",200,150))) +
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5)
p4

p<-ggpubr::ggarrange(p1,p2,p2a,p3,p4,nrow=5,ncol=1,heights=c(1.7,1.7,1.7,1.1,0.9),labels=c("a","b","c","d","e"))
ggsave(paste0(finalFigDir,"/",enhancerSet,"Enhancer_Promoter_contactRatio",fragLibs,"_maxK",maxKNN,"_fountVnonfount.pdf"),
       p,device="pdf",height=40, width=13, units="cm")



