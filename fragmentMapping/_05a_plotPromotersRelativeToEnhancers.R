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
# fountainsDir=paste0(projectDir,"/fountains")
# rnaSeqDir=paste0(projectDir,"/RNAseq_DGE")
# rnaTxSeqDir=paste0(projectDir,"/RNAseq_DTE")
# publicDataDir=paste0(projectDir,"/publicData")
# otherDataDir=paste0(projectDir,"/otherData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

source(paste0(projectDir,"/functions_plotting.R"))

# for large files from HiC fragment mapping
bigDataDir="/Volumes/external.data/MeisterLab/mdas/enhancer_fragment_mapping"
args=list()
## Parameters ----
if (length(args)==0) {
  minDistance=0
  maxDistance=30000
  enhancerSet="jaenes"
  #enhancerSet="daugherty"
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


#maxDistance=300000
maxKNN=5000
fragLibs="828v366"
# data from supl table S3 in paper.
totalTEVlib<-sum(258355543, 163424616)
totalCOH1lib<-sum(208298545, 130155406)
#totalSCC1lib<-sum(317973695, 250444362)
#totalSCC1COH1lib<-sum(172186520, 148507656)



## convert to table and analyse ------

fgi<-readRDS(paste0(bigDataDir,"/rds/04a__",enhancerSet,"Enh_Prom_",minDistance/1000,"-",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,"_prom",tssUpstream,"up",tssDownstream,"down.rds"))
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

# plots -----

### contact rank -----
pec_ss<-pec[pec$ctrlRank<=10,]
pec_ss$ctrlRank<-factor(pec_ss$ctrlRank,levels=c(1:10))
stats_df<-pec_ss %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(),ratioCoh1Ctrl=0)
stats_df
stat.test<-pec_ss  %>% wilcox_test(ratioCoh1Ctrl~ctrlRank,ref.group="1",
                                alternative="two.sided",
                                p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p1<-ggplot(pec_ss,aes(x=ctrlRank,y=ratioCoh1Ctrl)) +
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
pec_ss<-pec[pec$distanceRank<=10,]
pec_ss$distanceRank<-factor(pec_ss$distanceRank,levels=1:10)
stats_df<-pec_ss %>% dplyr::group_by(distanceRank) %>%
  dplyr::summarise(count=dplyr::n(),ratioCoh1Ctrl=0)
stats_df
stat.test<-pec_ss  %>% wilcox_test(ratioCoh1Ctrl~distanceRank,ref.group="1",
                                alternative="two.sided",
                                p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p2<-ggplot(pec_ss,aes(x=distanceRank,y=ratioCoh1Ctrl)) +
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
pec_ss<-pec[pec$ctrlRank<=10,]
pec_ss$ctrlRank<-factor(pec_ss$ctrlRank,levels=c(1:10))
stats_df<-pec_ss %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(),anchor2.LFC=-0.4)
stats_df
stat.test<-pec_ss %>% wilcox_test(anchor2.LFC~ctrlRank,ref.group="1",
                               alternative="two.sided",
                               p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p2a<-ggplot(pec_ss,aes(x=ctrlRank,y=anchor2.LFC)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_hline(yintercept=0,color="red")+
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  geom_signif(y_position=rep(c(0.33,0.38), nrow(stat.test))[1:nrow(stat.test)],
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0)+
  coord_cartesian(ylim=c(-0.4,0.4))  +
  ylab("RNAseq Log<sub>2</sub>FC") +
  xlab("Promoters ranked by contact probability in TEV control")
p2a

### LFC by distance -----
pec_ss<-pec[pec$distanceRank<=10,]
pec_ss$distanceRank<-factor(pec_ss$distanceRank,levels=c(1:10))
stats_df<-pec_ss %>% dplyr::group_by(distanceRank) %>% dplyr::summarise(count=dplyr::n(),anchor2.LFC=-0.4)
stats_df
stat.test<-pec_ss %>% wilcox_test(anchor2.LFC~distanceRank,ref.group="1",
                                  alternative="two.sided",
                                  p.adjust.method="fdr") %>%
  add_xy_position() %>% mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test

p2b<-ggplot(pec_ss,aes(x=distanceRank,y=anchor2.LFC)) +
  geom_boxplot(outlier.shape=NA,fill="grey90",notch=T) +
  geom_hline(yintercept=0,color="red")+
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
            size=2.5) +
  geom_signif(y_position=rep(c(0.33,0.38), nrow(stat.test))[1:nrow(stat.test)],
              annotations=stat.test$p.pretty,
              xmin=stat.test$group2,
              xmax=stat.test$group2,
              parse=T, size=0, textsize=2.5, tip_length=0)+
  coord_cartesian(ylim=c(-0.4,0.4))  +
  ylab("RNAseq Log<sub>2</sub>FC") +
  xlab("Promoters ranked by contact probability in TEV control")
p2b



pec_ss<-pec[pec$ctrlRank<=10,]
pec_ss$ctrlRank<-factor(pec_ss$ctrlRank,levels=1:10)
df1 = pec_ss[!duplicated(pec_ss$ctrlRank) ,]
p3<-ggplot(pec_ss,aes(x=distance/1000,group=ctrlRank)) + geom_histogram() +
  facet_wrap(.~ctrlRank,ncol=5) +
  xlab("Distance of promoter from enhancer (kb)")+
  ylab("Number of promoters ") +
  ggtitle("Contact probability rank in TEV control") +
  geom_text(data=df1,x=7.5,y=+Inf,aes(label=ctrlRank),
            vjust = "inward", hjust = "inward") +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
p3


pec_ss<-pec[pec$ctrlRank<=10,]
pec_ss$ctrlRank<-factor(pec_ss$ctrlRank,levels=1:10)
stats_df<-pec_ss %>% dplyr::group_by(ctrlRank) %>% dplyr::summarise(count=dplyr::n(), ctrlCounts=0)
stats_df
p4<-ggplot(pec_ss,aes(y=ctrlCounts,x=ctrlRank)) + geom_boxplot(outlier.shape=NA) +
  xlab("Promoter contact probability rank in TEV control")+
  ylab("Number of contacts ") +
  coord_cartesian(ylim=c(0,ifelse(enhancerSet=="daugherty",200,150))) +
  geom_text(data=stats_df,mapping=aes(label=count),color="blue",
                                            size=2.5)
p4

panelLabels<-if(enhancerSet=="daugherty"){c("a","c","e","g","i")} else {c("b","d","f","h","j")}
p<-ggpubr::ggarrange(p1,p2,p2a,p3,p4,nrow=5,ncol=1,heights=c(1,1,1,1.1,0.9),
                     labels=panelLabels)
ggsave(paste0(finalFigDir,"/",enhancerSet,"Enhancer_Promoter_contactRatio_",
            minDistance/1000,"-",maxDistance/1000,"kb_",fragLibs,"_maxK",maxKNN,
            "_prom",tssUpstream,"up",tssDownstream,"down.pdf"),
       p,device="pdf",height=30, width=13, units="cm")



