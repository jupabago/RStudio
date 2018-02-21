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

#This function takes the exact output from the matlab thresholded images
GetImageFiles<-function(folder, timepoints, positions){
  growthDf<-data.frame()
  fileList<-list.files(path = folder, full.names = TRUE)
  for (position in 0:positions){
    currentPosition<-fileList[grep(paste0("p",position,"_"),fileList)]
    for (timepoint in 0:timepoints){
      currentTimePoint<-currentPosition[grep(paste0("t",timepoint,"_"),currentPosition)]
      growthDf<-rbind(growthDf,NormRawPipeline(thresholdImage(currentTimePoint,"red", 512),thresholdImage(currentTimePoint,"green", 512),xysize=512, xydimSearch=25,sampleSize=1000,xyRealDim=1.24,pipeBinSize = 2,timePoint= timepoint, position = position))
    }}
  write.csv(growthDf, file = paste0(folder,"/PropOcData.csv"), row.names = FALSE)
  return(growthDf)
  }

PropOcCleanUp<-function(df){#This function removes irrelevant comparisons from the results dataframe
  df<-df[!(df$position==0 & df$source %in% c("finalRG","finalGG","finalGR")),]#remove first Sa mono
  df<-df[!(df$position==5 & df$source %in% c("finalRG","finalGG","finalGR")),]#remove second Sa mono
  df<-df[!(df$position==1 & df$source %in% c("finalRG","finalRR","finalGR")),]
  df<-df[!(df$position==6 & df$source %in% c("finalRG","finalRR","finalGR")),]
  df<-df[!(df$position==2 & df$source %in% c("finalRG","finalRR","finalGR")),]
  df<-df[!(df$position==7 & df$source %in% c("finalRG","finalRR","finalGR")),]
  return(df)
}
testCleanup<-PropOcCleanUp(testCleanup)
testTopHalf<-testCleanup[!(testCleanup$position>4),]
testBottomHalf<-testCleanup[!(testCleanup$position<5),]

popOcPA14_1<-GetImageFiles("/Volumes/Seagate Exp/GATECH images/11-14-17_timeSeries/images", timepoints = 10, positions =9)

GetGrowthCurve<-function(folder,timepoints, positions){
  growthDf<-data.frame()
  fileList<-list.files(path = folder, full.names = TRUE)
  for (position in 0:positions){
    currentPosition<-fileList[grep(paste0("p",position,"_"),fileList)]
    for (timepoint in 0:timepoints){
      currentTimePoint<-currentPosition[grep(paste0("t",timepoint,"_"),currentPosition)]
      growthDf<-rbind(growthDf,GetBiomass(thresholdImage(currentTimePoint,"red", 512),thresholdImage(currentTimePoint,"green", 512),xysize=512, xydimSearch=25,timePoint= timepoint, position))
    }
  }
  write.csv(growthDf, file = paste0(folder,"/GrowthCurve.csv", row.names = FALSE))
  return((growthDf))}

testGetData<-GetData(10)
testGetGrowthCurve<-GetGrowthCurve("/Volumes/Seagate Exp/GATECH images/11-8-17_timeSeries/images", 10, 9)
colnames(testGetGrowthCurve)<-c("Red", "Green", "TimePoint", "position")
ggplot()+
  geom_point(data = testGetGrowthCurve, aes(x=TimePoint, y = Green, color = "Pa", shape = factor(position)))+
  geom_point(data = testGetGrowthCurve, aes(x=TimePoint, y = Red, color = "Sa",shape = factor(position)))+
  scale_y_log10(name=expression(paste('Biomass' )), breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))

ggplot()+
  geom_point(data = testTopHalf, aes(x=distance, y = normPropOccup, color = source, shape = factor(position)))+
  #ggtitle(paste("Proportional occupancy r=",((((xydimSearch*2)+1)/2)*xyRealDim), sep = ''))+
  labs(title = "blah", y = "Proportional Occupancy", x = "Distance (um)")+
  scale_shape_discrete("condition", labels=c("Sa Mono", "PaWt Mono", "PaMut Mono", "PaWt Co", "PaMut Co"))+
  scale_color_discrete("relation", labels=c("Pa-Pa", "Sa-Pa", "Sa-Sa", "Pa-Sa"))+
  facet_grid(.~timePoint)
ggplot()+
  geom_point(data = testBottomHalf, aes(x=distance, y = normPropOccup, shape = source, color = factor(position)))+
  #ggtitle(paste("Proportional occupancy r=",((((xydimSearch*2)+1)/2)*xyRealDim), sep = ''))+
  labs(title = "blah", y = "Proportional Occupancy", x = "Distance (um)")+
  scale_color_discrete("condition", labels=c("Sa Mono", "PaWt Mono", "PaMut Mono", "PaWt Co", "PaMut Co"))+
  scale_shape_discrete("relation", labels=c("Pa-Pa", "Sa-Pa", "Sa-Sa", "Pa-Sa"))+
  facet_grid(.~timePoint)


