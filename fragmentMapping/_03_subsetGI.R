#!/usr/bin/env Rscript
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(InteractionSet)
library(plyranges)
library(GenomicInteractions)


args = commandArgs(trailingOnly=TRUE)

# for large files from HiC fragment mapping
bigDataDir="/mnt/external.data/MeisterLab/mdas/enhancer_fragment_mapping"
rnaSeqDir=paste0(bigDataDir,"/RNAseq_DGE")
publicDataDir=paste0(bigDataDir,"/publicData")


## Parameters ----
if (length(args)==0) {
  minDistance=0
  maxDistance=20000000
  tssUpstream=100
  tssDownstream=100
} else if (length(args)==4) {
  minDistance = as.numeric(args[1])
  maxDistance = as.numeric(args[2])
  tssUpstream = as.numeric(args[3])
  tssDownstream = as.numeric(args[4])
}

print(paste0("using command line args. minDistance:",minDistance/1000,"kb maxDistance:",maxDistance/1000,"kb"))
print(paste0("using ",tssUpstream,"bp upstream and ",tssDownstream,"bp downstream of tss"))

## Transcripts ----
gtfurl="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS19.canonical_geneset.gtf.gz"
gtffile<-gsub("\\.gz","",basename(gtfurl))

if(!file.exists(paste0(bigDataDir,"/",gtffile))){
  download.file(gtfurl,dest=basename(gtfurl))
  system(paste0("gunzip ",basename(gtfurl)))
}

if(!file.exists(paste0(bigDataDir,"/rds/03__transcriptPromoters_",tssUpstream,"upstream",tssDownstream,"downstream.rds"))){
  gtf<-import(paste0(bigDataDir,"/",gtffile))
  seqlevelsStyle(gtf)<-"UCSC"
  seqlevels(gtf)
  table(gtf$type)
  genes<-sort(gtf[gtf$type=="gene" & gtf$gene_biotype=="protein_coding"])
  txpt<-sort(gtf[gtf$type=="transcript" & gtf$gene_biotype=="protein_coding"])
  toBed<-txpt
  toBed$name<-txpt$transcript_id
  toBed$score<-1
  export(toBed,paste0(bigDataDir,"/bed/transcripts_WBPS19.bed"))
  txptPromoters<-promoters(txpt,upstream=tssUpstream,downstream=tssDownstream)
  saveRDS(txptPromoters,paste0(bigDataDir,"/rds/03__transcriptPromoters_",tssUpstream,"upstream",tssDownstream,"downstream.rds"))
}



## HiC fragments ----
fragments<-import(paste0(bigDataDir,"/bed/HicFrag_ce11.bed"))
strand(fragments)<-"*"



## fountains -----
#fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125c.RDS"))
#colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
#seqlevels(fountains)<-seqlevels(Celegans)
#fountains<-resize(fountains,width=6000,fix="center")



## RNAseq -------
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")
salmon<-readRDS(rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]



########################## pre processing #################

## join promoter fragments ------
if(!file.exists(paste0(bigDataDir,"/rds/03__joinedPromoterFragments_",tssUpstream,"up",tssDownstream,"down.rds"))){
  txptPromoters<-readRDS(paste0(bigDataDir,"/rds/03__transcriptPromoters_",tssUpstream,"upstream",tssDownstream,"downstream.rds"))
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

  promFrag<-data.frame(promFrag) %>% dplyr::group_by(seqnames,transcript_id) %>%
    dplyr::summarise(start=min(start),end=max(end),strand="*",
              gene_id=paste0(unique(gene_id),collapse=";"),
              fragNames=paste0(unique(name),collapse=";"),
              upVdown=paste0(unique(upVdown),collapse=";"),
              LFC=mean(LFC),
              padj=mean(padj,na.rm=T))%>% as_granges()
  table(promFrag$upVdown) # no mixed categories

  # eliminate duplicate fragment names by joining the relevant promoter names (same fragment
  # is promoter for two different genes)
  promFrag<-data.frame(promFrag) %>% dplyr::group_by(seqnames,fragNames) %>%
    summarise(start=min(start),end=max(end),strand="*",
              gene_id=paste0(unique(gene_id),collapse=";"),
              transcript_id=paste0(unique(transcript_id),collapse=";"),
              upVdown=paste0(unique(upVdown),collapse=";"),
              LFC=mean(LFC), padj=mean(padj)) %>% as_granges()
  table(promFrag$upVdown)
  # get rid of merged up/down and NS fragments
  promFrag$upVdown<-gsub(";NS|NS;","",promFrag$upVdown)
  #salmon[salmon$wormbaseID %in% c("WBGene00269391","WBGene00006975"),]
  promFrag$upVdown<-gsub("down;up","up",promFrag$upVdown) # only one and upreglated gene was more significant
  table(promFrag$upVdown)
  promFrag<-sort(promFrag)
  forBed<-promFrag
  forBed$name<-promFrag$fragNames
  export(forBed,paste0(bigDataDir,"/bed/joinedPromoterFragments_",tssUpstream,"up",tssDownstream,"down.bed"))
  saveRDS(promFrag,file=paste0(bigDataDir,"/rds/03__joinedPromoterFragments_",tssUpstream,"up",tssDownstream,"down.rds"))
}


## join enhancer fragments ----

for(enhancerSet in c("daugherty", "jaenes")){
  if(!file.exists(paste0(bigDataDir,"/rds/03__joined_",enhancerSet,"EnhancerFragments.rds"))){
    if(enhancerSet=="daugherty"){
      enhancers<-readRDS(paste0(publicDataDir,"/daugherty2017_L3enhancers_ce11.rds"))
      enhancers<-enhancers[enhancers$L3_chromHMMState=="L3_activeEnhancer",]
    }
    if(enhancerSet=="jaenes"){
      enhancers<-import(paste0(publicDataDir,"/Jaenes2018_L3activeEnhancers_ce11.bed"))
      enhancers$name<-paste0("peak",1:length(enhancers))
    }
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
    export(forBed,paste0(bigDataDir,"/bed/joined_",enhancerSet,"EnhancerFragments.bed"))
    saveRDS(enhFrag,paste0(bigDataDir,"/rds/03__joined_",enhancerSet,"EnhancerFragments.rds"))
  }
}




####################### subset fragments by enhancers and promoters
enhDaugherty<-readRDS(paste0(bigDataDir,"/rds/03__joined_daughertyEnhancerFragments.rds"))
enhJaenes<-readRDS(paste0(bigDataDir,"/rds/03__joined_jaenesEnhancerFragments.rds"))
promoters<-readRDS(paste0(bigDataDir,"/rds/03__joinedPromoterFragments_",tssUpstream,"up",tssDownstream,"down.rds"))

# combine all gr into one
allGR<-c(enhDaugherty,enhJaenes,promoters)
allGR<-unique(sort(allGR))



# subset ctrl
print("processing control gi")
ctrl<-readRDS(paste0(bigDataDir,"/rds/02__GI_366_fragment_pair_counts_",
              minDistance/1000,"-",maxDistance/1000,"kb.rds"))

# only keep interactions where both anchors overlap regions of interest
ctrlSubset<-ctrl[overlapsAny(ctrl,allGR,type="any",ignore.strand=T,use.region="first") &
      overlapsAny(ctrl,allGR,type="any",ignore.strand=T,use.region="second")]

saveRDS(ctrlSubset,paste0(bigDataDir,"/rds/03__GIss_366_fragment_pair_counts_",
                    minDistance/1000,"-",maxDistance/1000,"kb_enh_prom",tssUpstream,"up",tssDownstream,"down.rds"))


rm(ctrl)
rm(ctrlSubset)

# subset coh1
print("processing coh1 gi")
coh1<-readRDS(paste0(bigDataDir,"/rds/02__GI_828_fragment_pair_counts_",
                      minDistance/1000,"-",maxDistance/1000,"kb.rds"))

# only keep interactions where both anchors overlap regions of interest
coh1Subset<-coh1[overlapsAny(coh1,allGR,type="any",ignore.strand=T,use.region="first") &
                   overlapsAny(coh1,allGR,type="any",ignore.strand=T,use.region="second")]

saveRDS(coh1Subset,paste0(bigDataDir,"/rds/03__GIss_828_fragment_pair_counts_",
          minDistance/1000,"-",maxDistance/1000,"kb_enh_prom",tssUpstream,"up",tssDownstream,"down.rds"))

