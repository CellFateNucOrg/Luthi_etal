#! Rscript

library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(InteractionSet)
library(plyranges)
library(GenomicInteractions)

projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseq_DGE")
publicDataDir=paste0(projectDir,"/publicData")

bigDataDir="/mnt/external.data/MeisterLab/mdas/enhancer_fragment_mapping"


args = commandArgs(trailingOnly=TRUE)

## Parameters ----
if (length(args)==0) {
  minDistance=0
  maxDistance=20000000
  #enhancerSet="jaenes"
  enhancerSet="daugherty"
  tssUpstream=100
  tssDownstream=100
  print("using default parameters!")
  #stop("Please supply min and max distance for interaction", call.=FALSE)
} else if (length(args)==5) {
  minDistance = as.numeric(args[1])
  maxDistance = as.numeric(args[2])
  enhancerSet=args[3]
  tssUpstream=as.numeric(args[4])
  tssDownstream=as.numeric(args[5])
} else {
  print("need five ordered arguments: minDistance(bp)\n maxDistance(bp)\n enhancerSet('jaenes' or 'daugherty')\ntssUpstream(bp)\ntssDownstream(bp)\n")
}

print(paste0("using command line args. minDistance:",minDistance/1000,"kb maxDistance:",maxDistance/1000,"kb"))
print(paste0("using enhancer set: ",enhancerSet))
print(paste0("using promoters with ",tssUpstream,"bp upstream, and ",tssDownstream,"bp downstream of TSS"))


## Parameters ----
maxKNN=5000
fragLibs="828v366"
# data from supl table S3 in paper.
totalTEVlib<-sum(258355543, 163424616)
totalCOH1lib<-sum(208298545, 130155406)



## fountains -----
fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125c.RDS"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
seqlevels(fountains)<-seqlevels(Celegans)
fountains<-resize(fountains,width=6000,fix="center")



## RNAseq -------
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")
salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]
#sleuth<-readRDS("/Users/semple/Documents/MeisterLab/papers/Isiaka_etal/RNAseq_DTE/sleuth_coh1cs_DTE.RDS")



## functions ##########################
#gr1=promFrag
#gr2=enhFrag
#' function to find nearest k neighbours irrespective of orientation
#' and return GInteractions ofbject
findKNN<-function(gr1, gr2, minDistance, maxDistance, maxKNN){
  gr1w<-resize(gr1,width=maxDistance*2+1,fix="center")
  ol<-data.frame(findOverlaps(gr1w,gr2,ignore.strand=T))
  ol$distance<-GenomicRanges::distance(gr1[ol$queryHits],gr2[ol$subjectHits],ignore.strand=T)
  ol<-ol %>% dplyr::filter(distance>minDistance) %>%
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

## joined promoter fragments ------
promFrag<-readRDS(paste0(bigDataDir,"/rds/03__joinedPromoterFragments_",tssUpstream,"up",tssDownstream,"down.rds"))


## joined enhancer fragments ----
enhFrag<-readRDS(paste0(bigDataDir,"/rds/03__joined_",enhancerSet,"EnhancerFragments.rds"))


##  find all possible enhancer promoter pairs -----
fgi<-findKNN(enhFrag,promFrag, minDistance ,maxDistance,maxKNN)
table(fgi$distanceRank) # up to 783

# add pairwise distance of anchors
fgi$P_E_distanceMid<-pairdist(fgi)

#hist(fgi$P_E_distanceMid)

## add count data -----
if(!exists("ctrlgi")){
  ctrlgi<-readRDS(paste0(bigDataDir,"/rds/03__GIss_366_fragment_pair_counts_",
              minDistance/1000,"-",maxDistance/1000,"kb_enh_prom",tssUpstream,"up",tssDownstream,"down.rds"))
  ctrltotalFrag<-sum(ctrlgi$count)
}

ol<-data.frame(findOverlaps(fgi,ctrlgi))
ol$ctrlCounts<-ctrlgi$count[ol$subjectHits]
ol<-ol %>% dplyr::group_by(queryHits) %>% dplyr::mutate(sumCounts=sum(ctrlCounts))
fgi$ctrlCounts<-0
fgi$ctrlCounts[ol$queryHits]<-ol$sumCounts



if(!exists("coh1gi")){
  coh1gi<-readRDS(paste0(bigDataDir,"/rds/03__GIss_828_fragment_pair_counts_",
            minDistance/1000,"-",maxDistance/1000,"kb_enh_prom",tssUpstream,"up",tssDownstream,"down.rds"))
  coh1totalFrag<-sum(coh1gi$count)
}


ol<-data.frame(findOverlaps(fgi,coh1gi))
ol$coh1Counts<-coh1gi$count[ol$subjectHits]
ol<-ol %>% dplyr::group_by(queryHits) %>% dplyr::mutate(sumCounts=sum(coh1Counts))
fgi$coh1Counts<-0
fgi$coh1Counts[ol$queryHits]<-ol$sumCounts

# remove interactions wtih no counts in either dataset
sum(fgi$coh1Counts==0 & fgi$ctrlCounts==0)
fgi<-fgi[!(fgi$coh1Counts==0 & fgi$ctrlCounts==0)]


# label promoters with fountain location
fgi$fountVnonfountP<-"nonfount_P"
ol<-findOverlaps(anchors(fgi)$first,resize(fountains,width=60000,fix="center"))
fgi$fountVnonfountP[queryHits(ol)]<-"fountain_P"
ol<-findOverlaps(anchors(fgi)$first,fountains)
fgi$fountVnonfountP[queryHits(ol)]<-"tip_P"

# label enhancers with fountain location
fgi$fountVnonfountE<-"nonfount_E"
ol<-findOverlaps(anchors(fgi)$second,resize(fountains,width=60000,fix="center"))
fgi$fountVnonfountE[queryHits(ol)]<-"fountain_E"
ol<-findOverlaps(anchors(fgi)$second,fountains)
fgi$fountVnonfountE[queryHits(ol)]<-"tip_E"

# find distance and name of nearest fountain
fgi$anchor1fount<-fountains$fountainName[nearest(anchors(fgi)$first,fountains,ignore.strand=T)]
fgi$anchor2fount<-fountains$fountainName[nearest(anchors(fgi)$second,fountains,ignore.strand=T)]
fgi$anchor1fountDist<-mcols(distanceToNearest(anchors(fgi)$first,fountains,ignore.strand=T))$distance
fgi$anchor2fountDist<-mcols(distanceToNearest(anchors(fgi)$second,fountains,ignore.strand=T))$distance


# save
if(!file.exists(paste0(bigDataDir,"/rds/04a__",enhancerSet,"Enh_Prom_",minDistance/1000,"-",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,"_prom",tssUpstream,"up",tssDownstream,"down.rds"))){
  saveRDS(fgi,file=paste0(bigDataDir,"/rds/04a__",enhancerSet,"Enh_Prom_",minDistance/1000,"-",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,"_prom",tssUpstream,"up",tssDownstream,"down.rds"))
  export.bedpe(fgi, fn = paste0(bigDataDir,"/bedpe/",enhancerSet,"Enh_Prom_366_",minDistance/1000,"-",maxDistance/1000,"kb_maxK",maxKNN,"_prom",tssUpstream,"up",tssDownstream,"down.bedpe"), score = "ctrlCounts")
  export.bedpe(fgi, fn = paste0(bigDataDir,"/bedpe/",enhancerSet,"Enh_Prom_828_",minDistance/1000,"-",maxDistance/1000,"kb_maxK",maxKNN,"_prom",tssUpstream,"up",tssDownstream,"down.bedpe"), score = "coh1Counts")
}



