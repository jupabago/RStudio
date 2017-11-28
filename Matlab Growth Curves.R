##functions to analyze data from Matlab threshilding and overlapping aggregates
library(magrittr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)

parentFolder1<-"/Volumes/Seagate Exp/GATECH images/11-8-17_timeSeries/results"
parentFolder2<-"/Volumes/Seagate Exp/GATECH images/11-11-17_timeSeries/results"
parentFolder3<-"/Volumes/Seagate Exp/GATECH images/11-12-17_timeSeries/results"
parentFolder4<-"/Volumes/Seagate Exp/GATECH images/11-14-17_timeSeries/results"#first get the parent directory, this needs to be erased
fileList<-list.files(path = parentFolder, full.names = TRUE)
positions<-9
timepoints<-10

GetData<-function(folder,positions, timepoints){
  growthDf<-data.frame()
  for (position in 0:positions){
    currentPosition<-fileList[grep(paste0("p",position,"."),fileList)]
    for (timepoint in 0:timepoints){
      currentTimePoint<-currentPosition[grep(paste0("t",timepoint,"_"),currentPosition)]
      growthDf<-rbind(growthDf,GetStack(currentTimePoint, position, timepoint))
      }
  }
  growthDf<-growthDf[!(growthDf$position=="0" & growthDf$Bug=="Pa"),]
  growthDf<-growthDf[!(growthDf$position=="5" & growthDf$Bug=="Pa"),]
  growthDf<-growthDf[!(growthDf$position=="1" & growthDf$Bug=="Sa"),]
  growthDf<-growthDf[!(growthDf$position=="6" & growthDf$Bug=="Sa"),]
  growthDf<-growthDf[!(growthDf$position=="2" & growthDf$Bug=="Sa"),]
  growthDf<-growthDf[!(growthDf$position=="7" & growthDf$Bug=="Sa"),]
  return(growthDf)
  }

GetStack<-function(file, position, timepoint){
  datos <-read.csv(file, header = FALSE)
  colnames(datos)<-c('Total', 'Red', 'Green', 'Both', 'z_stack')
  datos%>%summarize(sum(Red))->redBiomass
  datos%>%summarize(sum(Green))->greenBiomass
  green<-c(greenBiomass, position, timepoint, "Pa")
  red<-c(redBiomass, position, timepoint,"Sa")
  output<-data.frame(green, stringsAsFactors=FALSE)
  output<-unname(output)
  output<-rbind(output,red)
  colnames(output)<-c('Biomass', 'position', 'timePoint', 'Bug')
  return(output)
  }

plotHist<-function(df){
  ggplot()+
    geom_point(data=df,aes(x=timePoint, y=Biomass, color = Bug), size=2)+
    scale_y_log10(name=expression(paste('Biomass' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    ggtitle("Growth Curves")+
    facet_wrap(~position, nrow = 2)
}


plotHist(Data1<-GetData(parentFolder1, positions, timepoints))
  
