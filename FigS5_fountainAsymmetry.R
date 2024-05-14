library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
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
source(paste0(projectDir,"/functions_fountainPlots.R"))


fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
colnames(mcols(fountains))[colnames(mcols(fountains))=="name"]<-"fountainName"
fountains$AsymmetryTercile<-paste0("Asym.Q",fountains$AsymmetryTercile)
fountains$AsymmetryQuintile<-paste0("Asym.Q",fountains$AsymmetryQuintile)



## 366 TPM - 10 kb zoom----
tpm366<-import(paste0(rnaSeqDir,"/tracks/PMW366_TPM_avr.bed"))
tpm366$length<-width(tpm366)
tpm366tss<-resize(tpm366,width=1,fix="start")
df<-getUpDownstreamFountainData(tpm366tss,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=10000,binSize=2000)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)


df$score<-NULL
df<-left_join(df,data.frame(tpm366)[,c(6,7)],by=c("name"))
df$log2TPMtev<-log2(df$score+1)
table(df$binnedDistance,df$AsymmetryQuintile)
table(df$binnedDistance1)


ignoreBins<-c(">-10","-0",">10","Asym.Q2","Asym.Q3","Asym.Q4")

selected=!df$binnedDistance1 %in% ignoreBins & !(df$AsymmetryQuintile %in% ignoreBins) & df$log2TPMtev >1
countLabels_a<-df[selected,]%>%
  dplyr::group_by(binnedDistance1,AsymmetryQuintile,.drop=T) %>%
  dplyr::summarise(count=dplyr::n())

# calculate stats
stat_df<-df[!(df$AsymmetryQuintile %in% ignoreBins) & df$log2TPMtev >0,] %>%
  rstatix::group_by(AsymmetryQuintile,binnedDistance,drop=T) %>% filter(binnedDistance==2) %>%
  rstatix::wilcox_test(log2TPMtev~binnedDistance1)  %>% filter(group2 %in% 2:10) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group1 %in% ignoreBins) %>%
  add_xy_position(x="binnedDistance1") %>%
  mutate(xmin=xmin-1,xmax=xmax-2,y.position[order(AsymmetryQuintile,y.position)])
stat_df

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p1<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2TPMtev),
               fill=rep(c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
                        length(unique(countLabels_a$AsymmetryQuintile))),
               outlier.color=NA) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=14,size=3,colour="blue",
            angle=0) +
  xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+1)")+
  ggtitle("log2 TPM in TEV control animals by fountain asymmetry") +
  facet_grid(AsymmetryQuintile~.) + coord_cartesian(ylim=c(0,15))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
                     hide.ns = F, color="purple",bracket.size=0.5,y=12)
p1
stat_df$p.adj





## coh-1 TPM - 10 kb zoom----
tpmCOH1<-import(paste0(rnaSeqDir,"/tracks/PMW828_TPM_avr.bed"))

tpmCOH1$length<-width(tpmCOH1)
tpmCOH1tss<-resize(tpmCOH1,width=1,fix="start")
df<-getUpDownstreamFountainData(tpmCOH1tss,fountains)

df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=10000,binSize=2000)
table(df$binnedDistance)

df<-flipDownstreamOrientation(df)
table(df$binnedDistance1)


df$score<-NULL
df<-left_join(df,data.frame(tpmCOH1)[,c(6,7)],by=c("name"))


table(df$binnedDistance,df$AsymmetryQuintile)

df$log2TPMcoh1<-log2(df$score+1)
ignoreBins<-c(">-10","-0",">10","Asym.Q2","Asym.Q3","Asym.Q4")


selected=!df$binnedDistance1 %in% ignoreBins & !(df$AsymmetryQuintile %in% ignoreBins) & df$log2TPMcoh1 >1
countLabels_a<-df[selected,]%>%
  dplyr::group_by(binnedDistance1,AsymmetryQuintile,.drop=T) %>%
  dplyr::summarise(count=dplyr::n())

# calculate stats
stat_df<-df[!(df$AsymmetryQuintile %in% ignoreBins) & df$log2TPMcoh1 >0,] %>%
  rstatix::group_by(AsymmetryQuintile,binnedDistance,drop=T) %>% filter(binnedDistance==2) %>%
  rstatix::wilcox_test(log2TPMcoh1~binnedDistance1)  %>% filter(group2 %in% 2:10) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj") %>%
  filter(!group1 %in% ignoreBins) %>%
  add_xy_position(x="binnedDistance1") %>%
  mutate(xmin=xmin-1,xmax=xmax-2,y.position[order(AsymmetryQuintile,y.position)])

bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)

p2<-ggplot(df[selected,]) +
  geom_boxplot(mapping=aes(x=binnedDistance1,y=log2TPMcoh1),
               fill=rep(c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
                        length(unique(countLabels_a$AsymmetryQuintile))),
               outlier.color=NA) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
  geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=14,size=3,colour="blue",
            angle=0) +
  xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+1)")+
  ggtitle("log2 TPM in COH-1 cleaved animals by fountain asymmetry") +
  facet_grid(AsymmetryQuintile~.) + coord_cartesian(ylim=c(0,15))+
  stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
                     hide.ns = F, color="purple",bracket.size=0.5,y=12)
p2
stat_df



######################################
## chromatin states ---------
######################################

states<-import.bed(paste0(publicDataDir,"/chromStates_L3_Evans2016_ce11.bed"))
#domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
#seqlevels(domains)<-seqlevels(states)

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")


tiles<-tileGenome(seqlengths(Celegans)[1:6],tilewidth=2000,cut.last.tile.in.chrom = T)
tiles<-getUpDownstreamFountainData(tiles,fountains)

tiles<-makeGRangesFromDataFrame(tiles,keep.extra.columns = T)

tiles$binnedDistance<-binByDistance1(tiles$distanceToFountain,maxDist=10000,binSize =2000,fullTile=T)
tiles<-flipDownstreamOrientation(tiles)
table(tiles$binnedDistance1,tiles$fountainLocation)
tiles$AsymmetryQuintile<-factor(tiles$AsymmetryQuintile,levels=paste0("Asym.Q",1:5))

listdf<-list()
for(l in levels(tiles$fountainLocation)){
  for(a in levels(tiles$AsymmetryQuintile)){
    for(d in levels(tiles$binnedDistance1)[c(2:6,8:13)]){
      df1<-getStateOLtable(gr=tiles[tiles$AsymmetryQuintile==a & tiles$binnedDistance1==d],states)
      df1$fountainLocation<-l
      df1$AsymmetryQuintile<-a
      df1$binnedDistance1<-d
      listdf[[paste0(l,"_",a,"_",d)]]<-df1
    }
  }
}


df<-do.call(rbind,listdf)
df$binnedDistance1<-factor(df$binnedDistance1,levels=levels(tiles$binnedDistance1)[c(2:6,8:13)])

df$XvA<-factor(df$XvA)

bin0<-ceiling(length(levels(df$binnedDistance1))/2)

boxdf<-df %>% filter(XvA=="A",AsymmetryQuintile %in% c("Asym.Q1","Asym.Q5")) %>%
  group_by(binnedDistance1,AsymmetryQuintile) %>%
  mutate(stateWidth=sum(stateWidth),stateFrequency=sum(stateFrequency)) %>%
  filter(binnedDistance1=="0") %>% dplyr::select(-AsymmetryQuintile,-state,-fountainLocation) %>% distinct()

p3<-ggplot(df[df$XvA=="A" & df$AsymmetryQuintile %in% c("Asym.Q1","Asym.Q5"),],
           aes(x=binnedDistance1,y=stateWidth/1e6,fill=state)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_manual(values=stateClrs) +
  xlab("Binned distance (kb)")+
  theme(legend.position = "none") +
  facet_grid(rows=vars(AsymmetryQuintile))+
  ggtitle(paste0("Chromatin state width in autosomal fountains")) +ylab("Bin coverage x Mb")
p3<-p3 + geom_bar(data=boxdf,aes(x=binnedDistance1, y=stateWidth/1e6,fill=NULL),stat="identity",
                  alpha=0,color="black",linewidth=0.5)

p3



p<-ggarrange(p1,p2,p3,nrow=3)
p<-annotate_figure(p, top = text_grob("Lüthi et al., Figure S5", size = 14))
ggsave(filename=paste0(finalFigDir,"/FigS5_fountainAsymmetry.pdf"),
       p, device="pdf",height=29,width=15, units="cm")






# ## coh-1 LFC - 10 kb zoom----
# lfc<-readRDS(paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds"))
# lfc<-GRanges(lfc[!is.na(lfc$chr),])
#
# lfc$length<-width(lfc)
# lfctss<-resize(lfc,width=1,fix="start")
# df<-getUpDownstreamFountainData(lfctss,fountains)
#
# df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=10000,binSize=2000)
# table(df$binnedDistance)
#
# df<-flipDownstreamOrientation(df)
# table(df$binnedDistance1)
#
#
# df$score<-NULL
#
# table(df$binnedDistance,df$AsymmetryQuintile)
#
# df$log2lfc<-df$log2FoldChange
# ignoreBins<-c(">-10","-0",">10","Asym.Q2","Asym.Q3","Asym.Q4")
#
#
# selected=!df$binnedDistance1 %in% ignoreBins & !(df$AsymmetryQuintile %in% ignoreBins)
# countLabels_a<-df[selected,]%>%
#   dplyr::group_by(binnedDistance1,AsymmetryQuintile,.drop=T) %>%
#   dplyr::summarise(count=dplyr::n())
#
# # calculate stats
# stat_df<-df[!(df$AsymmetryQuintile %in% ignoreBins) & df$log2FoldChange >0,] %>%
#   rstatix::group_by(AsymmetryQuintile,binnedDistance,drop=T) %>% filter(binnedDistance==2) %>%
#   rstatix::t_test(log2FoldChange~binnedDistance1)  %>% filter(group2 %in% 2:10) %>%
#   adjust_pvalue(method="BH") %>%
#   add_significance("p.adj") %>%
#   filter(!group1 %in% ignoreBins) %>%
#   add_xy_position(x="binnedDistance1") %>%
#   mutate(xmin=xmin-1,xmax=xmax-2,y.position[order(AsymmetryQuintile,y.position)])
#
# bin0<-ceiling(length(unique(countLabels_a$binnedDistance1))/2)
#
# p6<-ggplot(df[selected,]) +
#   geom_boxplot(mapping=aes(x=binnedDistance1,y=log2FoldChange),
#                fill=rep(c(rep("lightblue",bin0-1),"pink",rep("lightblue",bin0-1)),
#                         length(unique(countLabels_a$AsymmetryQuintile))),
#                outlier.color=NA) +
#   #geom_beeswarm(mapping=aes(x=binnedDistance1,y=log2FoldChange),size=0.5,alpha=0.3)+
#   theme(axis.text.x = element_text(angle = 0, vjust = 1),strip.text = element_text(size=10)) +
#   geom_text(data=countLabels_a,aes(x=binnedDistance1,label=count),y=0.45,size=3,colour="blue",
#             angle=0) +
#   xlab("Distance from fountain bins (kb)") + ylab("Log2(TPM+")+
#   ggtitle("log2 fold change in COH-1 cleaved animals by fountain asymmetry") +
#   facet_grid(AsymmetryQuintile~.) + coord_cartesian(ylim=c(-0.1,0.5))+
#   stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,
#                      hide.ns = F, color="purple",bracket.size=0.5,y=0.4) +
#   geom_hline(yintercept=0,linetype="dashed",colour="grey20")
# p6
# stat_df
#
