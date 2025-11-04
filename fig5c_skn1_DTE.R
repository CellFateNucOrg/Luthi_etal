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



ggpubr::ggarrange(pnull,p5,nrow=1,ncol=2,widths=c(0.9,0.2),labels=c(" ","h "))

p<-annotate_figure(p, top = text_grob("Isiaka et al., Figure 5", size = 14))
ggsave(paste0(finalFigDir,"/fig4_RNAseq.pdf"), p, device="pdf",
       width=19,height=19, units="cm")
