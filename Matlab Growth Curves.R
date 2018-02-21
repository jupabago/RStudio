##functions to analyze data from Matlab threshilding and overlapping aggregates
library(magrittr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)

GetData<-function(folder,positions, timepoints){
  fileList<-list.files(path = folder, full.names = TRUE)
  growthDf<-data.frame()
  for (position in 0:positions){
    currentPosition<-fileList[grep(paste0("p",position,"."),fileList)]
    for (timepoint in 0:timepoints){
      currentTimePoint<-currentPosition[grep(paste0("t",timepoint,"_"),currentPosition)]
      growthDf<-rbind(growthDf,GetStack(currentTimePoint, position, timepoint))
      }
  }
  #remove other color from monocultures
  #growthDf<-growthDf[!(growthDf$position=="0" & growthDf$Bug=="Pa"),]
  #growthDf<-growthDf[!(growthDf$position=="1" & growthDf$Bug=="Sa"),]
  #growthDf<-growthDf[!(growthDf$position=="2" & growthDf$Bug=="Sa"),]
  #growthDf<-growthDf[!(growthDf$position=="5" & growthDf$Bug=="Pa"),]
  #growthDf<-growthDf[!(growthDf$position=="6" & growthDf$Bug=="Sa"),]
  #growthDf<-growthDf[!(growthDf$position=="7" & growthDf$Bug=="Sa"),]
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

GetPositions<-function(df){
  SaMono<-df[(df$position=="0"| df$position=="5"),]
  PaWtMono<-df[(df$position=="1"| df$position=="6"),]
  PapqsLMono<-df[(df$position=="2"| df$position=="7"),]
  WtCo<-df[(df$position=="3"| df$position=="8"),]
  PqsLCo<-df[(df$position=="4"| df$position=="9"),]
  results<-gdata::combine(SaMono,PaWtMono, PapqsLMono,WtCo, PqsLCo)
  colnames(results)<-c('Biomass', 'position', 'timePoint', 'Bug','condition')
  return(results)
}

GetRatio<-function(df){
  SaMono<-df[(df$position=="0"| df$position=="5"),]
  PaWtMono<-df[(df$position=="1"| df$position=="6"),]
  PapqsLMono<-df[(df$position=="2"| df$position=="7"),]
  WtCo<-df[(df$position=="3"| df$position=="8"),]
  PqsLCo<-df[(df$position=="4"| df$position=="9"),]
  monoWtRatio<-PaWtMono$Biomass/SaMono$Biomass
  monoPqsLRatio<-PapqsLMono$Biomass/SaMono$Biomass
  coWtRatio<-WtCo[(WtCo$Bug=="Pa"),]$Biomass/WtCo[(WtCo$Bug=="Sa"),]$Biomass
  coPqsLRatio<-PqsLCo[(PqsLCo$Bug=="Pa"),]$Biomass/PqsLCo[(PqsLCo$Bug=="Sa"),]$Biomass
  ratiosDf<-data.frame(monoWtRatio,monoPqsLRatio,coWtRatio,coPqsLRatio, timepoint = sequence(11))
  results<-gather(ratiosDf, condition, timepoint)
  colnames(results)<-c("timepoint", "condition","ratio" )
  return(results)
}

plotGrowthCurve<-function(df, graphTitle, graphSubtitle = ''){
  ggplot()+
    geom_point(data=df,aes(x=timePoint, y=Biomass, color = Bug, shape = source), size=2)+
    scale_y_log10(name=expression(paste('Biomass' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    labs(title= graphTitle, subtitle = graphSubtitle)+
    facet_wrap(~condition, nrow = 2)
}

plotRatios<-function(df, colour, graphTitle, graphSubtitle = ''){
  if (colour=="source"){
    ggplot()+
      geom_point(data=df,aes(x=timepoint, y=ratio, color = source, shape = condition), size=2)+
      scale_y_log10(name=expression(paste('Biomass ratio (Pa/Sa)' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      labs(title= graphTitle, subtitle = graphSubtitle)
  }
  else {
    ggplot()+
      geom_point(data=df,aes(x=timepoint, y=ratio, color = condition, shape = source), size=2)+
      scale_y_log10(name=expression(paste('Biomass ratio (Pa/Sa)' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      labs(title= graphTitle, subtitle = graphSubtitle)
  }
}

plotMeans<-function(df, graphTitle, graphSubtitle = ''){
  ggplot(df, aes(x=timePoint, y = cellCount, color = condition, shape = Bug))+
    geom_errorbar(aes(ymin=cellCount-cellSD, ymax=cellCount+cellSD), width=.1) +
    geom_line() +
    geom_point()+
    labs(title= graphTitle, subtitle = graphSubtitle)+
    scale_y_log10(name=expression(paste('Biomass' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))
}