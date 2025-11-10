library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(rcompanion)
library(coin)
library(rstatix)
library(cowplot)
library(ggpubr)
library(S4Vectors)


options(tibble.width = Inf)
theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=10),
          axis.title.x=element_text(size=10))
)

projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseeq_DGE")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

source(paste0(projectDir,"/functions_plotting.R"))

#base_path<-"/Volumes/external.data/MeisterLab/pmeister/skn-1_deletions/240411_skn-1_deletions"
dataset="240411"
base_path=paste0(projectDir,"/skn1ImageAnalysis")

# ubs<-data.frame(ubs=c(" wt","ubs62","ubs71","ubs72","ubs73"),
#                 strain=c("1116","1182","1217",
#                          "1218","1219"),
#                 names=c(" wt","bec-1\u0394","nhr-46\u0394",
#                         "clec-178\u0394","1st intron\u0394"),
#                 order=c(1,4,5,2,3))
# ubs$combinedNames<-paste0(ubs$names,"\n(",ubs$ubs,")")

ubs<-data.frame(ubs=c("wt","ubs58","ubs59","ubs60","ubs62","ubs71","ubs72","ubs73",
                      "ubs72,ubs75","ubs72,ubs76"),
                strain=c("1116","1146","1163","1177","1182","1217",
                         "1218","1219","1241","1242"),
                names=c("wt","1st intron short\u0394","1st intron long\u0394",
                        "bec-1 short\u0394","bec-1\u0394","nhr-46\u0394",
                        "clec-178\u0394","1st intron\u0394",
                        "clec-178\u0394\n1st intron\u0394",
                        "clec-178\u0394\nbec-1\u0394"))
ubs$combinedNames<-paste0(ubs$names,"\n(",ubs$ubs,")")

selected<-c("wt","ubs71","ubs62","ubs73","ubs72","ubs59","ubs58","ubs60",
            "ubs72,ubs75","ubs72,ubs76")

idx<-match(selected,ubs$ubs)

ubs<-ubs %>% filter(ubs %in% selected)
ubs$combinedNames<-factor(ubs$combinedNames,levels=ubs$combinedNames[idx])



# scoring of nuclei

asi<-read.csv(list.files(paste0(base_path,"/imageStats/",dataset),pattern="^blobStats.*",full.names = T))

head(asi)

asifilt<-asi %>% dplyr::group_by(base_name) %>% dplyr::mutate(count=dplyr::n()) %>% dplyr::filter(count<=4)
#asifilt %>% dplyr::group_by(base_name) %>% dplyr::mutate(count=dplyr::n()) %>% filter(count>4) %>% group_by(base_name) %>%
#  dplyr::summarise(count=dplyr::n())
asifilt$strain<-factor(sapply(strsplit(asifilt$base_name,"_"),"[[",1))
asifilt$HSvNHS<-factor(sapply(strsplit(asifilt$base_name,"_"),"[[",2),levels=c("nHS","HS"))
levels(asifilt$HSvNHS)<-c("Control","Heatshock")
asifilt<-left_join(asifilt,ubs,by="strain")
asifilt$combinedNames<-factor(asifilt$combinedNames,levels=ubs$combinedNames)

sum.stat<- asifilt %>% dplyr::group_by(HSvNHS,combinedNames) %>%
  summarise(nuclei=dplyr::n(),avrIntensity=median(intensity_mean))
wtControl<-sum.stat %>% filter(HSvNHS=="Control", combinedNames=="wt\n(wt)")
wtHeatshock<-sum.stat %>% filter(HSvNHS=="Heatshock", combinedNames=="wt\n(wt)")
sum.stat$dividedByWtControl<-sum.stat$avrIntensity/wtControl$avrIntensity
sum.stat$dividedByWtHeatshock<-sum.stat$avrIntensity/wtHeatshock$avrIntensity
sum.stat$withHSvsWithout<-sum.stat$avrIntensity/sum.stat$avrIntensity[sum.stat$HSvNHS=="Control"]
sum.stat

stat.test <- asifilt %>% dplyr::group_by(HSvNHS) %>%
  rstatix::wilcox_test(intensity_mean~combinedNames,ref.group=levels(asifilt$combinedNames)[1],
                       p.adjust.method="fdr") %>%
  add_xy_position(x = "combinedNames")%>%
  mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj.signif!="ns")
stat.test$intensity_mean=1800
stat.test

stat.test2 <- asifilt %>% group_by(combinedNames) %>%
  rstatix::wilcox_test(intensity_mean~HSvNHS) %>%
  adjust_pvalue(method="fdr") %>%
  add_xy_position(x = "HSvNHS") %>%
  mutate(p.pretty = prettyExponents(p.adj,html=F)) %>%
  filter(p.adj<0.05)
stat.test2

facetlabels=c("no TEV induction", "cohesin<sup>COH-1</sup> cleavage")
names(facetlabels)=c("Control","Heatshock")


p5<-ggplot(asifilt,aes(x=combinedNames,y=intensity_mean)) +
  geom_boxplot(outlier.shape=NA,notch=T, fill="grey90") +
  facet_wrap(.~HSvNHS,labeller=as_labeller(facetlabels)) +
  geom_jitter(color="black",alpha=0.5,size=0.1,width=0.1) +
  theme(axis.title.x=element_blank())+
  ylab("Mean intensity (A.U.)") +
  coord_cartesian(ylim=c(0,2100))+
  geom_signif(data=stat.test,
              mapping=aes(x=group1,y=intensity_mean,group=HSvNHS),
              y_position=1800,
              annotations=stat.test$p.pretty,
              comparisons=zipup(stat.test$group1,stat.test$group2),
              parse=T, size=0.5, textsize=3, tip_length=0.01,colour="purple",
              manual=F)+
  geom_text(aes(label=nuclei),data=sum.stat,y=300,colour="black",size=3) +
  theme(panel.spacing = unit(1.3, "cm"),
        strip.text = ggtext::element_markdown(),
        axis.text.x=ggtext::element_markdown(angle=45,hjust=1)) +
  geom_hline(yintercept=sum.stat$avrIntensity[sum.stat$HSvNHS=="Control"& sum.stat$combinedNames=="wt\n(wt)"],
             color="green")+
  geom_hline(yintercept=sum.stat$avrIntensity[sum.stat$HSvNHS=="Heatshock"& sum.stat$combinedNames=="wt\n(wt)"],
             color="orange")
p5



# stat.test2$intensity_mean<- 100
# stat.test2$HSvNHS<-factor(rep("Heatshock",nrow(stat.test2)), levels=levels(asifilt$HSvNHS))
# p5<-p5+ geom_text(data=stat.test2,mapping=aes(x=combinedNames,y=intensity_mean,label=p.pretty),
#               size=3,vjust=1,parse=T,color="purple")

p5
p<-cowplot::plot_grid(p5,labels=c("e "))
p
ggsave(paste0(finalFigDir,"/Fig5e_skn1enhDelIntensity.pdf"), p,
       device=cairo_pdf,width=16,height=10,units="cm")


# compare variances
# https://www.datanovia.com/en/lessons/homogeneity-of-variance-test-in-r/
bartlett.test(intensity_mean ~ interaction(HSvNHS,combinedNames), data=asifilt)
fligner.test(intensity_mean ~ interaction(HSvNHS,combinedNames), data=asifilt)





