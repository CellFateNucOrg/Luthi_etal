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
library(tidyr)
library(purrr)
library(sleuth)


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


# skn-1 DTE ------

so<-readRDS(paste0(rnaTxSeqDir,"/sleuthObject.rds"))

counts1<-sleuth_to_matrix(so,"obs_norm","est_counts")
tmp<-counts1[grep("T19E7.2",rownames(counts1))[1:2],]
tmpdf<-data.frame(tmp)
tmpdf$transcript<-row.names(tmp)
tmpdf<-pivot_longer(tmpdf,cols=c(-transcript),names_to="sampleID",values_to="est_counts")

tmpdf$sample<-gsub("^X","",sapply(strsplit(tmpdf$sampleID,"_"),"[[",1))
tmpdf$replicate<-sapply(strsplit(tmpdf$sampleID,"_"),"[[",2)
tmpdf$lane<-sapply(strsplit(tmpdf$sampleID,"_"),"[[",3)

cnts366<-tmpdf %>% group_by(transcript,sample,replicate) %>%
  dplyr::summarise(average=mean(est_counts)) %>%
  filter(sample=="366")

cnts828<-tmpdf %>% group_by(transcript,sample,replicate) %>%
  dplyr::summarise(average=mean(est_counts)) %>%
  filter(sample=="828")

lfcpnts<-data.frame(b=log2(cnts828$average/cnts366$average),
                    label=rep(c("a","b"),each=3))


res<-data.frame(readRDS(paste0(rnaTxSeqDir,"/sleuth_coh1cs_DTE.RDS")))
skn1<-res$target_id[res$publicID=="skn-1"]
tmp<-res[res$publicID=="skn-1",]
tmp<-tmp[complete.cases(tmp),]
tmp$txt<-round(tmp$qval,2)
tmp$txt[tmp$txt>0.05]<-"ns"
tmp$txtPos<-0.24
tmp$label<-c("a","b")

p5<-ggplot(tmp,aes(x=label,y=b)) +
  geom_bar(stat="identity",fill="lightgrey",color="black") +
  xlab("*skn-1* isoform") + ylab(label="Log<sub>2</sub>FC") +
  geom_errorbar(aes(ymin=b-se_b,ymax=b+se_b),width=0.2) +
  geom_hline(yintercept=0) + coord_cartesian(ylim=c(-0.57,0.24)) +
  geom_text(aes(y=txtPos,label=txt),size=3,color="purple")+
  geom_jitter(data=lfcpnts,aes(x=label,y=b),width=0.1)
p5


ggsave(paste0(finalFigDir,"/fig5c_skn1_DTE.pdf"), p5, device="pdf",
       width=3,height=1.5, units="cm")
