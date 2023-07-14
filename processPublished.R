suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(BSgenome.Celegans.UCSC.ce11))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(kit))
suppressPackageStartupMessages(library(Organism.dplyr))


projectDir="."
rnaSeqDir=paste0(projectDir,"/RNAseq_DGE")
publicDataDir=paste0(projectDir,"/publicData")

genomeDir=paste0("/Users/semple/Documents/MeisterLab/GenomeVer/",genomeVer)
genomeVer="WS285"


source(paste0(rnaSeqDir,"/functions.R"))
source(paste0(rnaSeqDir,"/DESeq2_functions.R"))

if(!dir.exists(publicDataDir)) {
  dir.create(publicDataDir)
}


ce11seqinfo<-seqinfo(Celegans)
#metadata<-readRDS(paste0(rnaSeqDir,"/wbGeneGR_WS275.rds"))


ce10toCe11url<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
ce10toCe11<-"ce10Toce11.over.chain"
download.file(ce10toCe11url,paste0(projectDir,"/publicData/",ce10toCe11,".gz"))
system(paste0("gunzip ",projectDir,"/publicData/",ce10toCe11,".gz"))
chainCe10toCe11<-import.chain(paste0(projectDir,"/publicData/",ce10toCe11))

# gene metadata
metadata<-getMetadataGR(genomeDir,genomeVer,projectDir,format="ce11")

#########################-
# Enhancer datasets -----
#########################-

#########################-
## Daugherty ATAC enhancers -----
#########################-

#Genome Res. 2017 Dec;27(12):2096-2107.
# doi: 10.1101/gr.226233.117. Epub 2017 Nov 15.
# Chromatin accessibility dynamics reveal novel functional enhancers in C. elegans
# Aaron C Daugherty 1, Robin W Yeo 1, Jason D Buenrostro 1, William J Greenleaf 1 2, Anshul Kundaje 1 3, Anne Brunet 1 4
# https://pubmed.ncbi.nlm.nih.gov/29141961/
# For all analyses, the ce10/WS220 version of the C. elegans genome (Rosenbloom et al. 2015)was used.

if(!file.exists(paste0(projectDir,"/publicData/daugherty2017_enhancers_ce11.rds"))){
  daughertyURL1<-"https://genome.cshlp.org/content/suppl/2017/11/15/gr.226233.117.DC1/Supplemental_Table_S3.xlsx"
  #daughertyURL2<-"https://genome.cshlp.org/content/suppl/2017/11/15/gr.226233.117.DC1/Supplemental_Table_S10.xlsx"

  download.file(daughertyURL1,destfile=paste0(projectDir,"/publicData/",basename(daughertyURL1)))
  #download.file(daughertyURL2,destfile=paste0(projectDir,"/publicData/",basename(daughertyURL2)))

  daughertytab3<-read_excel(paste0(projectDir,"/publicData/",basename(daughertyURL1)))
  #daughertytab10<-read_excel(paste0(projectDir,"/publicData/",basename(daughertyURL2)))

  table(daughertytab3$L3_chromHMMState)

  daughertyEnhGR<-GRanges(daughertytab3)
  daughertyEnhGR
  #liftover to ce11
  daughertyEnhGR_ce11<-unlist(liftOver(daughertyEnhGR,chain=chainCe10toCe11))
  saveRDS(daughertyEnhGR_ce11,paste0(projectDir,"/publicData/Daugherty2017_SuplTbl_S3_ce11.rds"))

  n<-GenomicRanges::nearest(daughertyEnhGR_ce11,metadata)
  daughertyEnhGR_ce11$nearestGene<-metadata$wormbaseID[n]
  daughertyEnhGR_ce11$distanceToNearest<-mcols(distanceToNearest(daughertyEnhGR_ce11,metadata))$distance
  daughertyEnhGR_ce11

  daughertyEnhL3<-daughertyEnhGR_ce11[daughertyEnhGR_ce11$L3_chromHMMState %in%
                                        c("L3_activeEnhancer",
                                          "L3_repressedEnhancer",
                                          "L3_H3K27me3Repressed"),]

  hist(daughertyEnhL3$distanceToNearest,breaks=100)
  saveRDS(daughertyEnhL3,file=paste0(projectDir,"/publicData/daugherty2017_L3enhancers_ce11.rds"))


  daughertyEnhEE<-daughertyEnhGR_ce11[daughertyEnhGR_ce11$EE_chromHMMState %in%
                                        c("EE_activeEnhancer",
                                          "EE_repressedEnhancer",
                                          "EE_H3K27me3Repressed"),]
  hist(daughertyEnhEE$distanceToNearest,breaks=100)
  saveRDS(daughertyEnhEE,file=paste0(projectDir,"/publicData/daugherty2017_EEenhancers_ce11.rds"))

  daughertyEnhYA<-daughertyEnhGR_ce11[daughertyEnhGR_ce11$YA_chromHMMState %in%
                                        c("YA_activeEnhancer",
                                          "YA_repressedEnhancer",
                                          "YA_H3K27me3Repressed"),]
  hist(daughertyEnhYA$distanceToNearest,breaks=100)
  saveRDS(daughertyEnhYA,file=paste0(projectDir,"/publicData/daugherty2017_YAenhancers_ce11.rds"))

  file.remove(paste0(projectDir,"/publicData/",basename(daughertyURL1)))
}



#########################-
## Jaenes ATAC enhancers -----
#########################-
# Chromatin accessibility dynamics across C. elegans development and ageing
# Jürgen Jänes, Yan Dong, Michael Schoof, Jacques Serizay, Alex Appert, Chiara Cerrato, Carson Woodbury, Ron Chen, Carolina Gemma, Ni Huang, Djem Kissiov, Przemyslaw Stempor, Annette Steward, Eva Zeiser, Sascha Sauer, Julie Ahringer
# Oct 26, 2018
# https://doi.org/10.7554/eLife.37344


#' Get indeces of GRanges that have biggest overlap
#'
#' Takes in two GRanges objects and finds the indeces of the second that have the
#' longest overlap with the first, and returns and indexing object for both.
#' @param gr1 First genomic ranges
#' @param gr2 Second genomic ranges
#' @return Table of indeces of gr1 and gr2 with the gr2 ranges that have longest overlap with gr1
#' @export
whichBiggestOverlap<-function(gr1,gr2){
  gr1$uniqRangeIDgr1<-paste0(seqnames(gr1),":",start(gr1),"-",end(gr1))
  gr2$uniqRangeIDgr2<-paste0(seqnames(gr2),":",start(gr2),"-",end(gr2))
  length(gr1)
  grol<-join_overlap_intersect(gr1,gr2)
  grol$olWidth<-width(grol)
  groll<-grol %>% group_by(uniqRangeIDgr1) %>% summarise(numOverlap=n(),overlapSize=width, overlapIDgr1=uniqRangeIDgr1,
                                                         overlapIDgr2=uniqRangeIDgr2)
  groll$orderedIDs<-groll$overlapIDgr2[sapply(groll$overlapSize,order)]
  groll$topID<-mapply(function(idList,last) idList[[last]],groll$orderedIDs,groll$numOverlap)
  groll$overlapIDgr1<-sapply(groll$overlapIDgr1,"[[",1)

  idx1<-match(groll$overlapIDgr1,gr1$uniqRangeIDgr1)
  idx2<-match(groll$topID,gr2$uniqRangeIDgr2)
  df<-data.frame(GR1=idx1,GR2=idx2)
  return(df)
}

#' Add Daugherty et al. ChromHMM states to a genomic ranges
#'
#' Add Daugherty et al.(Genome Res, 2017) ChromHMM states (input as a GRanges object)
#' to another gr. If multiple states overlap, the function will take
#' the longest, unless it is "Low signal", in which case it will take
#' the second longest more informative state.
#' @param gr1 Genomic ranges object for which you want to find chromHMM state
#' @param chromHMMgr ChromHMM GRanges object from Daugherty et al. specific to EEmb, L3 or YAd
#' @param stage String indicating which stage the chromHMM is taken from.
#' @return returns GRanges (gr1) with additional columns for "numStateOverlap",
#' "overlapSize","overlapState" ,"topState". "topState" Contains the chosen state
#' @export
addBiggestOverlapChromHMM<-function(gr1,chromHMMgr,stage){
  gr1$uniqRangeID<-paste0(seqnames(gr1),":",start(gr1),"-",end(gr1))
  ji<-join_overlap_intersect(gr1,chromHMMgr)
  ji$olWidth<-width(ji)
  jii<-ji %>% group_by(uniqRangeID) %>% summarise(numStateOverlap=n(),overlapSize=width, overlapState=state)
  jii$orderedStates<-jii$overlapState[sapply(jii$overlapSize,order)]
  jii$topState<-mapply(function(stateList,last) stateList[[last]],jii$orderedStates,jii$numStateOverlap)
  # replace "Low signal" by next best overlap if possible
  jii$topState[jii$topState=="Low signal" & jii$numStateOverlap>1]<-mapply(function(stateList,last) stateList[[last-1]],
                                                                           jii$orderedStates[jii$topState=="Low signal" & jii$numStateOverlap>1],
                                                                           jii$numStateOverlap[jii$topState=="Low signal" & jii$numStateOverlap>1])
  jii$orderedStates<-NULL
  jiigr<-GRanges(jii$uniqRangeID)
  colnames(jii)<-c("uniqRangeID", paste0(colnames(jii)[-1],"_",stage,"_chromHMM"))
  mcols(jiigr)<-jii
  gr11<-GRanges(left_join(data.frame(gr1),data.frame(jii),by=c("uniqRangeID")))
  return(gr11)
}




if(!file.exists(paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11.rds"))){
  jaenesURL1<-"https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzczNDQvZWxpZmUtMzczNDQtZmlnMi1kYXRhMS12Mi50eHQ-/elife-37344-fig2-data1-v2.txt?_hash=jnh09dk%2F9t%2BIseamB5NWBCgtLxFmYQ%2BPJIOMmxucAww%3D"
  jaenesFile1<-"elife-37344-fig2-data1-v2.txt"

  download.file(jaenesURL1,jaenesFile1)
  jaenes1<-read.delim(jaenesFile1,sep="\t",header=T)
  dim(jaenes1)
  head(jaenes1)
  table(jaenes1$annot)
  jaenes1<-jaenes1[jaenes1$annot=="putative_enhancer",]
  jaenes1gr<-GRanges(seqnames=jaenes1$chrom_ce10,ranges=IRanges(start=jaenes1$start_ce10,end=jaenes1$end_ce10),
                     strand="*")
  jaenes1gr$uniqRangeID<-paste0(seqnames(jaenes1gr),":",start(jaenes1gr),"-",end(jaenes1gr))
  # read in atac peak data
  jaenesURL3<-"https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzczNDQvZWxpZmUtMzczNDQtZmlnMS1kYXRhMS12Mi50eHQ-/elife-37344-fig1-data1-v2.txt?_hash=S6tOfhfttAOIosQvxZIyFwFKgJ5yfYfpnXuvXIwyv%2BQ%3D"
  jaenesFile3<-"elife-37344-fig1-data1-v2.txt"
  download.file(jaenesURL3,jaenesFile3)
  jaenes3<-read.delim(jaenesFile3,sep="\t",header=T)
  dim(jaenes3)
  jaenes3gr<-GRanges(jaenes3)
  # merge datasets
  idx<-whichBiggestOverlap(jaenes1gr,jaenes3gr)
  jaenes<-cbind(jaenes1[idx$GR1,c(4:9,16:17,20:21)],jaenes3[idx$GR2,c(4:9)])
  stageNames<-names(jaenes[grep("atac_",names(jaenes))])
  stageNamesShort<-gsub("_height$","",gsub("atac_wt_","",stageNames))
  # find stage at which peak is highest, and some metrics of its specificty
  jaenes$maxStage<-stageNamesShort[apply(jaenes[,stageNames],1,which.max)]
  jaenes$maxStagePeak<-apply(jaenes[,stageNames],1,max)
  jaenes$ratioToSecond<-round(jaenes$maxStagePeak/sapply(apply(jaenes[,stageNames],1,topn,n=2L,index=F,simplify=F),min),2)
  jaenes$peakSD<-round(apply(jaenes[,stageNames],1,sd),3)
  jaenesgr<-GRanges(seqnames=jaenes$chrom_ce11,
                    ranges=IRanges(start=jaenes$start_ce11,end=jaenes$end_ce11),
                    strand="*")
  mcols(jaenesgr)<-jaenes[5:20]
  jaenesgr<-sort(sortSeqlevels(jaenesgr))
  saveRDS(jaenesgr,paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11_stages.rds"))
  # # export bedGraph files with atac signal for each stage
  # for(s in 1:length(stageNames)){
  #   forBG<-GRanges(seqnames=jaenes$chrom_ce11,
  #                  ranges=IRanges(start=jaenes$start_ce11,end=jaenes$end_ce11),
  #                  strand="*")
  #   forBG$score<-jaenes[,stageNames[s]]
  #
  #   export(forBG,paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11_",
  #                       stageNamesShort[s],".bedGraph"))
  # }
  # add chromHMM data from Daugherty2017
  jaenesgr<-readRDS(paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11_stages.rds"))
  # L3 stage
  chromHMM<-readRDS("/Users/semple/Documents/MeisterLab/Datasets/Daugherty2017_GenomeRes_GSE89608/chromHMM_L3_Daugherty2017_ce11.rds")
  jaenesgr<-addBiggestOverlapChromHMM(jaenesgr,chromHMM,stage="L3")
  #EEmb
  chromHMMemb<-readRDS("/Users/semple/Documents/MeisterLab/Datasets/Daugherty2017_GenomeRes_GSE89608/chromHMM_EEmb_Daugherty2017_ce11.rds")
  jaenesgr<-addBiggestOverlapChromHMM(jaenesgr,chromHMMemb,stage="EEmb")
  saveRDS(jaenesgr,paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11_stages_chromHMM.rds"))
  # make bed files for QC
  forBed<-jaenesgr
  mcols(forBed)<-NULL
  forBed$name<-jaenesgr$topState_L3_chromHMM
  export(forBed,paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11_stages_L3chromHMM.bed"))
  forBed$name<-jaenesgr$topState_EEmb_chromHMM
  export(forBed,paste0(projectDir,"/publicData/Jaenes2018_enhancers_ce11_stages_EEmbchromHMM.bed"))
  # clean up
  file.remove(jaenesFile1)
  file.remove(jaenesFile3)
}



#####################-
## Chromatin states Evans et al. (2016) (Ahinger lab)------
#####################-
# https://www.pnas.org/content/pnas/113/45/E7020.full.pdf

# Dataset S1. Coordinates of EE and L3 chromatin states
# File of chromosome, start position, end position, and state number. Coordinates are in WS220 and follow BED conventions (start positions are in zero-based coordinates and end positions in one-based coordinates).

if( !file.exists(paste0(projectDir,"/publicData/chromDomains_L3_Evans2016_ce11.bed"))){
  chromStatesURL<-"https://www.pnas.org/doi/suppl/10.1073/pnas.1608162113/suppl_file/pnas.1608162113.sd01.xlsx"
  download.file(chromStatesURL,paste0(projectDir,"/publicData/",basename(chromStatesURL)))
  ce10toCe11url<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
  ce10toCe11<-"ce10Toce11.over.chain"
  download.file(ce10toCe11url,paste0(projectDir,"/publicData/",ce10toCe11,".gz"))
  system(paste0("gunzip ",projectDir,"/publicData/",ce10toCe11,".gz"))
  #file.remove(paste0(projectDir,"/publicData/",ce10toCe11,".gz"))

  chrAstates<-readxl::read_excel(paste0(projectDir,"/publicData/pnas.1608162113.sd01.xlsx"),sheet="L3 autosome states",col_names=c("chr","start","end","state"))
  chrXstates<-readxl::read_excel(paste0(projectDir,"/publicData/pnas.1608162113.sd01.xlsx"),sheet="L3 chr X states",col_names=c("chr","start","end","state"))


  chrStates<-rbind(chrAstates,chrXstates)
  chrStatesGR<-GRanges(paste0("chr",chrStates$chr,":",
                              (chrStates$start+1),"-",
                              chrStates$end))
  chrStatesGR$score<-c(chrAstates$state,chrXstates$state)
  chainCe10toCe11<-import.chain(paste0(projectDir,"/publicData/",ce10toCe11))
  chrStatesGR_ce11<-unlist(liftOver(chrStatesGR,chain=chainCe10toCe11))

  seqinfo(chrStatesGR_ce11)<-ce11seqinfo
  export(chrStatesGR_ce11,
         con=paste0(projectDir,"/publicData/chromStates_L3_Evans2016_ce11.bed"),
         format="bed")
  file.remove(paste0(projectDir,"/publicData/",basename(chromStatesURL)))
}

