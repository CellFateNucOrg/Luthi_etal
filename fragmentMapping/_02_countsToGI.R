#!/usr/bin/env Rscript
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(InteractionSet)
library(plyranges)
library(GenomicInteractions)
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)

# for large files from HiC fragment mapping
bigDataDir="/mnt/external.data/MeisterLab/mdas/enhancer_fragment_mapping"


## Parameters ----
if (length(args)==0) {
  minDistance=0
  maxDistance=20000000
  #stop("Please supply min and max distance for interaction", call.=FALSE)
} else if (length(args)==2) {
  minDistance = as.numeric(args[1]) 
  maxDistance = as.numeric(args[2]) 
  print(paste0("using command line args. minDistance:",minDistance/1000,"kb maxDistance:",maxDistance/1000,"kb")) 
}


## HiC fragments ----
fragments<-import(paste0(bigDataDir,"/bed/HicFrag_ce11.bed"))
strand(fragments)<-"*"



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


###########################


## get control counts ----
if(!file.exists(paste0(bigDataDir,"/rds/02__GI_366_fragment_pair_counts_",
                       minDistance/1000,"-",maxDistance/1000,"kb.rds"))){
  pairCount1<-read.delim(paste0(bigDataDir,"/fragCounts/366_3_fragment_pair_counts_",
                                minDistance/1000,"-", maxDistance/1000,"kb.txt"),
                         sep=" ",header=F)
  counts1<-sum(pairCount1$V2)
  counts1
  gi1<-makeCisGI(pairCount1,fragments,maxDistance)
  rm(pairCount1)
  pairCount2<-read.delim(paste0(bigDataDir,"/fragCounts/366_4_fragment_pair_counts_",
                                minDistance/1000,"-", maxDistance/1000,"kb.txt"),
                         sep=" ",header=F)  #163424616
  counts2<-sum(pairCount2$V2)
  counts2
  gi2<-makeCisGI(pairCount2,fragments,maxDistance)
  rm(pairCount2)
  ctrlgi<-combineGIcounts(gi1,gi2)
  rm(gi1,gi2)
  #forBedpe<-ctrlgi
  #forBedpe$name<-paste0(ctrlgi$anchor1.name,"|",ctrlgi$anchor2.name)
  #forBedpe$score<-ctrlgi$count
  #export(forBedpe,paste0("GI_366_fragment_pair_counts_",minDistance/1000,
  #"-",maxDistance/1000,"kb.bedpe"))
  saveRDS(ctrlgi,paste0(bigDataDir,"/rds/02__GI_366_fragment_pair_counts_",
                        minDistance/1000,"-",maxDistance/1000,"kb.rds"))
}


## get coh-1 counts ----
if(!file.exists(paste0(bigDataDir,"/rds/02__GI_828_fragment_pair_counts_",
                       minDistance/1000,"-",maxDistance/1000,"kb.rds"))){
  pairCount1<-read.delim(paste0(bigDataDir,"/fragCounts/828_1_fragment_pair_counts_",
                                minDistance/1000,"-",maxDistance/1000,"kb.txt"),
                         sep=" ",header=F)
  gi1<-makeCisGI(pairCount1,fragments,maxDistance)
  rm(pairCount1)
  pairCount2<-read.delim(paste0(bigDataDir,"/fragCounts/828_2_fragment_pair_counts_",
                                minDistance/1000,"-",maxDistance/1000,"kb.txt"),
                         sep=" ",header=F)
  gi2<-makeCisGI(pairCount2,fragments,maxDistance)
  rm(pairCount2)
  coh1gi<-combineGIcounts(gi1,gi2)
  rm(gi1,gi2)
  #forBedpe<-coh1gi
  #forBedpe$name<-paste0(coh1gi$anchor1.name,"|",coh1gi$anchor2.name)
  #forBedpe$score<-coh1gi$count
  #export(forBedpe,paste0("GI_828_fragment_pair_counts_",maxDistance/1000,"kb.bedpe"))
  saveRDS(coh1gi,paste0(bigDataDir,"/rds/02__GI_828_fragment_pair_counts_",
                        minDistance/1000,"-",maxDistance/1000,"kb.rds"))
}

