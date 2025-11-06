library(readxl)
library(tidyr)
library(lattice)
library(nlme)
library(splines)

#excel tables: rows represent individual worms, columns larval stage (1-4) or molt (hatch, M1-4)
#ldH: larval stage duration in hours
#muv: relative growth rate in 1/h
#volAtEcdysis: volume of individual worms at ecdysis: hatch, M1-M5 (for hatch only the worms hatched in the imaging AFTER HS have the volume)
#validworms: worms that hatched at least 2.5 h before the HS

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

workPath="/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains"
experiments<-data.frame(date=c(20230210,20230505),path=paste0(workPath,"/growthAssay/",c("20230210__Results_of_COH-1_TEV_heat-shock_experiment","20230505__Results_of_COH-1_TEV_heat-shock_experiment")))

# read in data
df=NULL
for( e in 1:nrow(experiments)){
  for(strain in 1:2){
    valid<-read_excel(paste0(experiments$path[e],"/validworms",strain,".csv"),
                      col_names=F,
                      range = cellranger::cell_limits(ul =c(1L, 1L)))
    colnames(valid)<-c("validWorms")
    ecdys<-read_excel(paste0(experiments$path[e],"/volAtEcdysis",strain,".csv"),
                    col_names=c("H",paste0("M",1:4)),
                    range = cellranger::cell_limits(ul =c(1L, 1L)))
    tmp<-data.frame(date=experiments$date[e],strain=strain,validWorms=valid[,1])
    tmp<-cbind(tmp,ecdys)
    if(is.null(df)){
      df<-tmp
    } else {
      df<-rbind(df,tmp)
    }
  }
}



df$treatment<-factor(c("coh-1cs","control")[df$strain], levels=c("control","coh-1cs"),
                     labels=c("TEV control","COH-1\ncleavage"))
df$worm<-paste0("worm",1:nrow(df))
df$H<-NULL
dfkeep<-df

# keep valid worms
df<-df[df$validWorms,]

head(df)
tail(df)
dim(df)

#### larval volume at ecdysis------
vol<-df
table(vol$date,vol$treatment)
vol<-pivot_longer(vol,cols=grep("M.{1}",colnames(vol)),names_to="stage",values_to="volume")
vol$stage<-factor(vol$stage)
vol$worm<-factor(vol$worm)
vol$date<-factor(vol$date)
print(vol,width=Inf)
dim(vol)



bxp<-ggplot(data=vol,aes(x=stage,y=volume, color=treatment)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(alpha=0.5,size=1,position = position_jitterdodge()) +
  facet_wrap(~as.numeric(date)) +
  scale_y_continuous(trans='log10') +
  ylab("Volume at moult (nanoliters)") +
  labs(color="Animals")

stat.test<-vol %>% dplyr::group_by(date,stage) %>% rstatix::wilcox_test(volume~treatment) %>%
  adjust_pvalue() %>%  add_significance("p.adj")
print(stat.test,width=Inf)

stat.test <- stat.test %>% add_x_position(x="stage")
stat.test <- stat.test %>% add_y_position(y.trans=log10)
#print(stat.test,width=Inf)
stat.test$p.format <- p_format(
  stat.test$p.adj, accuracy = 0.001,
  leading.zero = T
)
stat.test$custom.label <- ifelse(stat.test$p.adj <= 0.05, stat.test$p.format, "ns")


p<-bxp +
  stat_pvalue_manual(stat.test, label = "custom.label",tip.length=0)
p
ggsave(paste0(finalFigDir,"/FigS13c_growthAssay.pdf"),p, width=16,height=10,units="cm")

