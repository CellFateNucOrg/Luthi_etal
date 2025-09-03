library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)

ubs<-data.frame(ubs=c("wt","ubs58","ubs59","ubs60","ubs62","ubs71","ubs72","ubs73"),
                strain=c("1116","1146","1163","1177","1182","1217",
                         "1218","1219"))

projectDir="."
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

#base_path<-"Z:/MeisterLab"
base_path<-"/Volumes/external.data/MeisterLab"


dataset<-"240411"
scorefilename1<-paste0(base_path,"/jsemple/skn1_deletions/",dataset,"_ImageScore_v2.xlsx")
scorefilename2<-paste0(base_path,"/jsemple/skn1_deletions/",dataset,"_ImageScore_PM.xlsx")
datafilename<-paste0(base_path, "/pmeister/skn-1_deletions/",dataset,"_skn-1_deletions/imageStats/allImageStats_0.5_99.95_169_702.csv")

dataset<-"240410"
scorefilename1<-paste0(base_path,"/jsemple/skn1_deletions/",dataset,"_ImageScore_v1.xlsx")
scorefilename2<-paste0(base_path,"/jsemple/skn1_deletions/",dataset,"_ImageScore_PM.xlsx")
datafilename<-paste0(base_path, "/pmeister/skn-1_deletions/",dataset,"_skn-1_deletions/imageStats/allImageStats_0.5_99.95_142_434.csv")


sheetNum=2
df<-read.csv(datafilename)
scores1<-read_excel(scorefilename1, sheet=sheetNum)
scores2<-read_excel(scorefilename2, sheet=sheetNum)
colnames(scores2)<-c("blind_key","score_PM","twoWorms")
allScores<-left_join(scores1,scores2)
allScores$twoWorms<-ifelse(is.na(allScores$twoWorms),F,T)
allScores<-rbind(allScores, allScores[allScores$twoWorms==T,])

df1<-left_join(allScores,df,by=join_by("blind_key"))

df1$meanScore<-round(rowMeans(cbind(df1$score,df1$score_PM)),0)

df1$strain<-factor(sapply(strsplit(df1$base_name,"_"),"[[",1))
df1$HSvNHS<-factor(sapply(strsplit(df1$base_name,"_"),"[[",2),levels=c("nHS","HS"))
df1$meanScore<-factor(df1$meanScore,levels=c(1:4))
df1<-left_join(df1,ubs,by="strain")
df1$ubs<-factor(df1$ubs,levels=ubs$ubs)

df2<-df1 %>% group_by(ubs, HSvNHS) %>% dplyr::summarise(count=n())
df2$ubs<-factor(df2$ubs)
df2$meanScore<-NA

p<-ggplot(df1,aes(x=ubs,group=meanScore,fill=meanScore)) +
  geom_bar(position = position_fill(reverse=T),color="black",linewidth=0.1) +
  facet_wrap(.~HSvNHS) +
  scale_fill_grey(start=1,end=0,na.translate=F) +
  geom_text(data=df2,aes(x=ubs,label=count),y=1.025,  color="blue",
            size=3)
p

ggsave(paste0(finalFigDir,"/supplFig_skn1Del_AvrScores_",dataset,".pdf"), p,
       device="pdf",width=14,height=8,units="cm")

df1 %>% group_by(HSvNHS, strain) %>% mutate(totalCount=n()) %>%
  group_by(HSvNHS, strain, meanScore) %>% summarise(count=n(),fraction=count/unique(totalCount))

asi<-read.csv(paste0(base_path,"/pmeister/skn-1_deletions/240411_skn-1_deletions/imageStats/blobStats_0.5_99.95_169_702.csv"))
dim(asi)
asi<-read.csv(paste0(base_path,"/pmeister/skn-1_deletions/240410_skn-1_deletions/imageStats/blobStats_0.5_99.95_142_434.csv"))

asi<-read.csv(paste0(base_path,"/pmeister/skn-1_deletions/240501_skn-1_deletions/imageStats/blobStats_0.5_99.95_171_816.csv"))
head(asi)
asifilt<-asi %>% group_by(base_name) %>% mutate(count=n()) %>% filter(count<=2)

asifilt$strain<-factor(sapply(strsplit(asifilt$base_name,"_"),"[[",1))
asifilt$HSvNHS<-factor(sapply(strsplit(asifilt$base_name,"_"),"[[",2),levels=c("nHS","HS"))

ggplot(asifilt,aes(x=strain,y=intensity_mean)) +
  geom_boxplot(outlier.shape=NA,notch=T, fill="grey80") +
  facet_wrap(.~HSvNHS) + geom_beeswarm(color="grey20",alpha=0.5)

ggplot(asifilt,aes(x=strain,y=area)) +
  geom_boxplot(outlier.shape=NA,notch=T, fill="grey80") +
  facet_wrap(.~HSvNHS) + geom_beeswarm(color="grey20",alpha=0.5)

