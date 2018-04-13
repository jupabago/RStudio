library(rdist)
library(pixmap)
library(rtiff)
library(reshape2)
library(gplots)
library(dplyr)
library(gdata)
library(ggplot2)
library(cowplot)

source('/Users/jupabago/testRstudioGithub/PropOcFunctionsPipeline.R')
source('/Users/jupabago/testRstudioGithub/Matlab Growth Curves.R')

#----------------------------------GROWTH CURVES--------------------------
dataPA14_exp_1<-GetData("/Volumes/Seagate Exp/GATECH images/2-8-18_timeSeries/results", positions = 7, timepoints = 43)
dataPA14_exp_2<-GetData("/Volumes/Seagate Exp/GATECH images/2-9-18_timeSeries/results", positions = 7, timepoints = 43)
#dataWLMpao1_3<-GetData("/Volumes/Seagate Exp/GATECH images/11-2-17_timeSeries/results", positions = 9, timepoints = 9) need to rewrite "GetData" function

GetPositions1<-function(df){#this function exists in "MatlabGrowthCurves.R" but a bit different
  WtCo<-df[(df$position=="0"| df$position=="2"| df$position=="4"| df$position=="6"),]
  PqsLCo<-df[(df$position=="1"| df$position=="3"| df$position=="5"| df$position=="7"),]
  results<-gdata::combine(WtCo, PqsLCo)
  colnames(results)<-c('Biomass', 'position', 'timePoint', 'Bug','condition')
  return(results)
}

plotMeans1<-function(df, graphTitle, graphSubtitle = ''){#this function exists in "MatlabGrowthCurves.R" but a bit different
  ggplot()+
    geom_errorbar(aes(ymin=cellCount-cellSD, ymax=cellCount+cellSD), width=.1) +
    geom_line(data = df, aes(x=timePoint, y = cellCount, color = Bug)) +
    geom_point(data = df, aes(x=timePoint, y = cellCount, color = Bug, shape = condition))+
    labs(title= graphTitle, subtitle = graphSubtitle)+
    scale_y_log10(name=expression(paste('Biomass' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))
}

GetGrowthCurve<-function(df){
allData<-GetPositions1(df)
allData%>% group_by(Bug, condition, timePoint) %>% summarise(cellCount = mean(Biomass), cellSD = sd(Biomass)) %>% as.data.frame() -> allDataMean
allDataMean%>%filter(timePoint<22)->allDataMean
#allDataMean$timePoint<-factor(allDataMean$timePoint)
#levels(allDataMean$timePoint)<-c("hour 2","hour 2.5","hour 3","hour 3.5","hour 4","hour 4.5","hour 5","hour 5.5","hour 6","hour 6.5","hour 7","hour 7.5","hour 8","hour 8.5","hour 9","hour 9.5","hour 10","hour 10.5","hour 11","hour 11.5","hour 12","hour 12.5")
#plotMeans1(allDataMean, "Pa Sa co-culture Growth Curve 1")
ggplot()+
  geom_errorbar(aes(x=timePoint,ymin=cellCount-cellSD, ymax=cellCount+cellSD,color = Bug),data = allDataMean,width=.1) +
  geom_line(data = allDataMean, aes(x=timePoint, y = cellCount, color = Bug, linetype = condition)) +
  geom_point(data = allDataMean,size = 3, aes(x=timePoint, y = cellCount, color = Bug, shape = condition))+
  labs(title= "Pa Sa co-culture Growth Curve 1")+
  scale_y_log10(name=expression(paste('Biomass' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))
}
GetGrowthCurve(dataPA14_exp_1)
GetGrowthCurve(dataPA14_exp_2)


#---------------------------------SPATIAL ORGANIZATION-------------------
#
#the following two lines collect the spatial organization data that was output of the "GetImageFiles" function

#This function assigns the correct condition to the different positions
GetPositions2<-function(df){#this function exists in "MatlabGrowthCurves.R" but a bit different
  WtCo<-df[(df$position=="0"| df$position=="2"| df$position=="4"| df$position=="6"),]
  PqsLCo<-df[(df$position=="1"| df$position=="3"| df$position=="5"| df$position=="7"),]
  results<-gdata::combine(WtCo, PqsLCo)
  return(results)
}
#Read csv files generated from script...
pa14_exp_1<-read.csv("~/Documents/Whiteley/PROJECTS/pqsL SCFM2/Prop oc time series data/Matlab thresholding/PA14_exp_1_PropOcData.csv", header = TRUE)
pa14_exp_2<-read.csv("~/Documents/Whiteley/PROJECTS/pqsL SCFM2/Prop oc time series data/Matlab thresholding/PA14_exp_2_PropOcData.csv", header = TRUE)

pa14_exp_1_clean<-GetPositions2(pa14_exp_1)
pa14_exp_2_clean<-GetPositions2(pa14_exp_2)

graphStatsPropOc<-function(df){
df%>%select(distance, normPropOccup, timePoint, position, source, source.1)%>%filter(source %in% c("finalGR","finalRG"))%>%filter(timePoint<22)%>%group_by(distance,timePoint,source, source.1)%>%summarise(avgPropOc=mean(normPropOccup),SDPropOc=sd(normPropOccup) )%>% as.data.frame()->statsPropOc
statsPropOc$timePoint<-factor(statsPropOc$timePoint)
levels(statsPropOc$timePoint)<-c("hour 2","hour 2.5","hour 3","hour 3.5","hour 4","hour 4.5","hour 5","hour 5.5","hour 6","hour 6.5","hour 7","hour 7.5","hour 8","hour 8.5","hour 9","hour 9.5","hour 10","hour 10.5","hour 11","hour 11.5","hour 12","hour 12.5")
ggplot()+
  geom_point(data = statsPropOc, aes(x=distance, y = avgPropOc, color = source, shape= source.1), size = 3)+
  geom_line(data = statsPropOc, aes(x=distance, y = avgPropOc, color = source, shape= source.1))+
  geom_errorbar(aes(x=distance, ymin=avgPropOc-SDPropOc, ymax=avgPropOc+SDPropOc,color=source), data = statsPropOc,width=.1) +
  labs(title = paste0("Proportional Occupancy over time 1"), y = "Proportional Occupancy", x = "Distance (um)", shape= "mutant", color="Bug")+
  facet_wrap(~timePoint, nrow = 2)+
  scale_colour_discrete(name="Focus->Surrounding",breaks=c("finalGR", "finalRG"),labels=c("Pa->Sa", "Sa->Pa"))+
  scale_shape_discrete(name = "Condition", breaks=c("WtCo", "PqsLCo"), labels=c("Sa+PA14wt", "Sa+PA14-pqsL"))+
  scale_size_continuous(guide = FALSE)+
  theme(legend.position="bottom")
}
graphStatsPropOc(pa14_exp_1_clean)
graphStatsPropOc(pa14_exp_2_clean)



#Empty Space functions
GraphEmptySpaceTimeCourse1<-function(df,title){#this function exists in "MatlabGrowthCurves.R" but a bit different
  df%>%filter(timePoint<22)->df
  df$timePoint<-factor(df$timePoint)
  levels(df$timePoint)<-c("hour 2","hour 2.5","hour 3","hour 3.5","hour 4","hour 4.5","hour 5","hour 5.5","hour 6","hour 6.5","hour 7","hour 7.5","hour 8","hour 8.5","hour 9","hour 9.5","hour 10","hour 10.5","hour 11","hour 11.5","hour 12","hour 12.5")
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(0,2,4,6))%>%filter(source %in% c("finalRR","finalRG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->rwCSat
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(0,2,4,6))%>%filter(source %in% c("finalGR","finalGG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->gwCSat
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(1,3,5,7))%>%filter(source %in% c("finalRR","finalRG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->rmCSat
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(1,3,5,7))%>%filter(source %in% c("finalGR","finalGG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->gmCSat
  ggplot()+
    geom_line(data = gwCSat, aes(x=distance, y = AvgEmpty, color = "Pa"))+
    geom_point(data = gwCSat,size =3, aes(x=distance, y = AvgEmpty, color = "Pa", shape="wt"))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Pa"), data = gwCSat,width=.1) +
    
    geom_line(data = gmCSat, aes(x=distance, y = AvgEmpty, color = "Pa"))+
    geom_point(data = gmCSat,size =3, aes(x=distance, y = AvgEmpty, color = "Pa", shape="-pqsL"))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Pa"), data = gmCSat,width=.1) +
    
    geom_line(data = rwCSat, aes(x=distance, y = AvgEmpty, color = "Sa"))+
    geom_point(data = rwCSat,size =3, aes(x=distance, y = AvgEmpty, color = "Sa", shape="wt"))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Sa"), data = rwCSat,width=.1) +
    
    geom_line(data = rmCSat, aes(x=distance, y = AvgEmpty, color = "Sa"))+
    geom_point(data = rmCSat,size =3, aes(x=distance, y = AvgEmpty, color = "Sa", shape="-pqsL"))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Sa"), data = rmCSat,width=.1) +
    
    labs(title = paste0("Empty space over time",title), y = "Proportional Occupancy", x = "Distance (um)", color="condition")+
    facet_wrap(~timePoint, nrow = 2)}

GraphEmptySpaceTimeCoursePa<-function(df,title){#this function exists in "MatlabGrowthCurves.R" but a bit different
  df%>%filter(timePoint<22)->df
  df$timePoint<-factor(df$timePoint)
  levels(df$timePoint)<-c("hour 2","hour 2.5","hour 3","hour 3.5","hour 4","hour 4.5","hour 5","hour 5.5","hour 6","hour 6.5","hour 7","hour 7.5","hour 8","hour 8.5","hour 9","hour 9.5","hour 10","hour 10.5","hour 11","hour 11.5","hour 12","hour 12.5")
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(0,2,4,6))%>%filter(source %in% c("finalGR","finalGG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame()->gwCSat
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(1,3,5,7))%>%filter(source %in% c("finalGR","finalGG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->gmCSat
  ggplot()+
    geom_line(data = gwCSat, aes(x=distance, y = AvgEmpty, color = "Pa wt co "))+
    geom_point(data = gwCSat,size=3, aes(x=distance, y = AvgEmpty, color = "Pa wt co "))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Pa wt co "), data = gwCSat,width=.1) +
    geom_line(data = gmCSat, aes(x=distance, y = AvgEmpty, color = "Pa mut co "))+
    geom_point(data = gmCSat,size=3, aes(x=distance, y = AvgEmpty, color = "Pa mut co "))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Pa mut co "), data = gmCSat,width=.1) +
    labs(title = paste0("Empty space over time around Pa",title), y = "Proportional Occupancy", x = "Distance (um)", color="condition")+
    facet_wrap(~timePoint, nrow = 2)}

GraphEmptySpaceTimeCourseSa<-function(df,title){#this function exists in "MatlabGrowthCurves.R" but a bit different
  df%>%filter(timePoint<22)->df
  df$timePoint<-factor(df$timePoint)
  levels(df$timePoint)<-c("hour 2","hour 2.5","hour 3","hour 3.5","hour 4","hour 4.5","hour 5","hour 5.5","hour 6","hour 6.5","hour 7","hour 7.5","hour 8","hour 8.5","hour 9","hour 9.5","hour 10","hour 10.5","hour 11","hour 11.5","hour 12","hour 12.5")
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(0,2,4,6))%>%filter(source %in% c("finalRR","finalRG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->rwCSat
  df%>%select(distance, normPropOccup, timePoint, position, source)%>%filter(position %in% c(1,3,5,7))%>%filter(source %in% c("finalRR","finalRG"))%>%group_by(distance,timePoint, position)%>%summarise(total=sum(normPropOccup))%>%mutate(emptySpace = 1-total)%>% group_by(distance, timePoint) %>% summarise(AvgEmpty = mean(emptySpace), SDEmpty = sd(emptySpace)) %>% as.data.frame() ->rmCSat
  ggplot()+
    geom_point(data = rwCSat, size = 3,aes(x=distance, y = AvgEmpty, color = "Sa wt co "))+
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Sa wt co "), data = rwCSat,width=.1) +
    geom_errorbar(aes(x=distance, ymin=AvgEmpty-SDEmpty, ymax=AvgEmpty+SDEmpty,color = "Sa mut co "), data = rmCSat,width=.1) +
    geom_line(data = rwCSat, aes(x=distance, y = AvgEmpty, color = "Sa wt co ")) +
    geom_point(data = rmCSat,size = 3, aes(x=distance, y = AvgEmpty, color = "Sa mut co "))+
    geom_line(data = rmCSat, aes(x=distance, y = AvgEmpty, color = "Sa mut co ")) +
    labs(title = paste0("Empty space over time around Sa",title), y = "Proportional Occupancy", x = "Distance (um)", color="condition")+
    facet_wrap(~timePoint, nrow = 2)}

GraphEmptySpaceTimeCourse1(pa14_exp_1_clean, " rep1")
GraphEmptySpaceTimeCourseSa(pa14_exp_1_clean, " rep1")
GraphEmptySpaceTimeCourseSa(pa14_exp_2_clean, " rep2")
GraphEmptySpaceTimeCoursePa(pa14_exp_1_clean, " rep1")
GraphEmptySpaceTimeCoursePa(pa14_exp_2_clean, " rep2")





