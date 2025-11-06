library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(rcompanion)
library(coin)
library(rstatix)
library(cowplot)

theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=10),
          axis.title.x=element_text(size=10))
)

projectDir="."
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}


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

alleles<-ubs

scores=avrScores

makeScorePlot<-function(datafile,scores,alleles,sheetNum=2){
  df<-datafile
  df<-left_join(df,scores,by=join_by("blind_key"))

  df$strain<-factor(sapply(strsplit(df$base_name,"_"),"[[",1))
  df<-left_join(df,alleles,by=join_by("strain"))

  worm2<-df[!is.na(df$score2),]
  worm2$score<-worm2$score2
  df<-rbind(df,worm2)

  df$HSvNHS<-factor(sapply(strsplit(df$base_name,"_"),"[[",2),levels=c("nHS","HS"))
  levels(df$HSvNHS)<-c("Control","Heatshock")
  df$score<-factor(df$score,levels=c("c",1:4))

  df<-df[!is.na(df$score),]
  df1<-df %>% group_by(combinedNames, HSvNHS) %>% dplyr::summarise(count=n())
  df1$combinedNames<-factor(df1$combinedNames,levels=alleles$combinedNames)
  df1$score<-NA

  p<-ggplot(df,aes(x=combinedNames,group=score,fill=score)) +
    geom_bar(position = position_fill(reverse=T),color="black",linewidth=0.1) +
    facet_wrap(.~HSvNHS) +
    scale_fill_grey(start=1,end=0,na.translate=F) +
    coord_cartesian(ylim=c(0,1.08))+
    geom_text(data=df1,aes(x=combinedNames,label=count),y=1.05,  color="blue",
              size=3) +
    theme(axis.text.x = element_text(angle=0,hjust=0.5,face="italic"),
          axis.title.x=element_blank())
  return(list(p, df))
}

scorers=c("JS","PM")
dataset<-"240411"
# dataset<-"240410"
# dataset<-"240501"
# dataset<-"240502"
# dataset<-"240515"
# dataset<-"240516"

base_path<-paste0("/Volumes/external.data/MeisterLab/pmeister/skn-1_deletions/",dataset,"_skn-1_deletions/")
datafile<-read.csv(list.files(paste0(base_path, "/imageStats/"), pattern="^allImageStats.*\\.csv$",full.names=T))


# combine scores
score1<-read_excel(paste0(base_path,"/scoring/",dataset,"_ImageScore_",scorers[1],".xlsx"),sheet=2)
score2<-read_excel(paste0(base_path,"/scoring/",dataset,"_ImageScore_",scorers[2],".xlsx"),sheet=2)


combineScores<-function(score1,score2){
  scoreavr<-score1[,1:3]
  scoreavr$score<-round(rowMeans(cbind(as.numeric(score1$score),as.numeric(score2$score)),na.rm=T),0)
  scoreavr$score2<-round(rowMeans(cbind(as.numeric(score1$score2),as.numeric(score2$score2)),na.rm=T),0)
  worm2<-scoreavr[!is.na(scoreavr$score2),]
  worm2$score<-worm2$score
  scoreavr<-rbind(scoreavr,worm2)
  scoreavr$score2<-NA
  scoreavr<-scoreavr[!is.na(scoreavr$score),]
  return(scoreavr)
}

avrScores<-combineScores(score1,score2)


p3<-makeScorePlot(datafile,avrScores,ubs)
p3[[1]]

results<-p3[[2]]%>% filter(HSvNHS=="Heatshock",score!="c") %>% group_by(combinedNames,score) %>%
  summarise(count=n())
results$score<-droplevels(results$score)

if(dataset=="240411"){
  ftab<-xtabs(count ~ combinedNames + score ,data=results)
  #prop.table(ftab)
  # Extended Cochran–Armitage test ----
  # for ordered nominal variable vs non ordered nominal variable
  ridx<-which(rownames(ftab) %in% levels(ubs$combinedNames)[1:5])
  coin::chisq_test(ftab[ridx,],scores=list("score"=c(1:4)))

  PT = pairwiseOrdinalIndependence(ftab[ridx,], compare="row")

  pvaldf<-data.frame(label=round(PT[1:4,"p.adjust"],2),
                     HSvNHS = factor("Heatshock",levels=c("Control","Heatshock")),
                     group1 =sapply(strsplit(PT[1:4,"Comparison"],":"),"[[",1),
                     group2=sapply(strsplit(PT[1:4,"Comparison"],": "),"[[",2))

  pvaldf$group2
  p3[[1]]<-p3[[1]]+annotate(geom="text",x=2,y=0.8,label=pvaldf$label[1],colour="white",size=3)
  p3[[1]]<-p3[[1]]+annotate(geom="text",x=3,y=0.8,label=pvaldf$label[2],colour="white",size=3)
  p3[[1]]<-p3[[1]]+annotate(geom="text",x=4,y=0.8,label=pvaldf$label[3],colour="white",size=3)
  p3[[1]]<-p3[[1]]+annotate(geom="text",x=5,y=0.8,label=pvaldf$label[4],colour="white",size=3)
}
p3[[1]]<-p3[[1]]+ylab("Fraction of worms")
p3

p<-cowplot::plot_grid(p3[[1]],labels=c("b "))
p
ggsave(paste0(finalFigDir,"/FigS12b_skn1enhDelScores.pdf"), p,
       device=cairo_pdf,width=19,height=10,units="cm")



