##----testa ----
library(rdist)
library(pixmap)
library(rtiff)
library(reshape2)
library(gplots)
library(dplyr)
library(gdata)
library(ggplot2)
library(cowplot)

#### 1. THE SETUP
rm(list = ls())

getFiles<-function(path){
  setwd(path)  
  parent <- getwd()
  files <- dir(parent, recursive = T)
  tif <- files[grep("tif", files)]
  fileList<-list(c1=tif[grep("Aa", tif)],c2=tif[grep("Sg", tif)])
}
thresholdImage<-function(filesList,channel,side){
  c1 <- array(0, c(side, side, length(filesList)))
  for(i in 1:length(filesList)){
    print(paste("loading ch1...", i, "/", length(filesList), filesList[i]))
    c1[,,i]<-pixmap::getChannels(readTiff(filesList[i]), colors = channel)
  }
  c1_t <- array(0, c(side, side, length(filesList)))
  for(i in 1:dim(c1)[3]){
    print(paste("thresholding ch1...", i, "/", length(filesList), filesList[i]))
    if(mean(c1[,,i])==0){next}#the code breaks if you try to apply "autothreshold"" on an empty image, so this prevents it
    th <- autoThreshold(c1[,,i], mean(c1[,,i]))[2]
    c1_t[,,i][c1[,,i] > th] <- 1 
  }
  return(c1_t)
}
CorrectSampleVolume<-function(thresholdedImage, xytotal, ztotal, xydim, zdim){#this limits the boundaries of voxel sampling to prevent going out of bounds
  correctedImage<- thresholdedImage[  thresholdedImage$x >= xydim & thresholdedImage$x <= (xytotal-xydim) &
                                        thresholdedImage$y >= xydim & thresholdedImage$y <= (xytotal-xydim) &
                                        thresholdedImage$z >= zdim & thresholdedImage$z <= (ztotal-zdim),]
  return(correctedImage)
}
IdVoxel<-function(thresholdedImage){#this makes a df with the coordinates from all voxel objects
  dfVoxelCoords <- data.frame(which(thresholdedImage == 1, T))
  colnames(dfVoxelCoords) <- c("x", "y", "z")
  return(dfVoxelCoords)
}
SampleVoxels<-function(voxelsCoords, sampleSize){#randomly samples non-empty voxels from a thresholded image
  sampledVoxels <- sample(1:dim(voxelsCoords)[1], size = sampleSize)
  sampleVoxelsCoords <- voxelsCoords[sampledVoxels,]
  return(sampleVoxelsCoords)
}

#Only one of the following 2 functions may be used, DistanceRadius has one extra line to filter out distances out of the radius
DistanceBox<-function(focus,thresholdedImage, xydim, zdim, xyscale, zscale){
  #takes a voxel and collects all the surrounding voxels within the specified box
  #input of xydim and zdim is in voxels
  outbox<-thresholdedImage[(focus$x-xydim):(focus$x+xydim), (focus$y-xydim):(focus$y+xydim), (focus$z-zdim):(focus$z+zdim)]
  outboxCoords<-(IdVoxel(outbox))#find non-empty voxels in box
  outboxCoords<-outboxCoords-1#move to the origin
  outboxCoords$x<-outboxCoords$x*xyscale#convert pixels to real lengths
  outboxCoords$y<-outboxCoords$y*xyscale
  outboxCoords$z<-outboxCoords$z*zscale
  distances<-as.data.frame(t(cdist(t(as.matrix(c(xydim*xyscale, xydim*xyscale, zdim*zscale))), as.matrix(outboxCoords))))
  return(distances)
}
DistanceRadius<-function(focus,thresholdedImage, xydim, zdim, xyscale, zscale){
  #takes a voxel and collects all the surrounding voxels within the specified sphere
  #input of xydim and zdim is in voxels
  #maximum distance allowed is capped at xydim*xyscale
  outbox<-thresholdedImage[(focus$x-xydim):(focus$x+xydim), (focus$y-xydim):(focus$y+xydim), (focus$z-zdim):(focus$z+zdim)]
  outboxCoords<-(IdVoxel(outbox))#find non-empty voxels in box
  outboxCoords<-outboxCoords-1#move to the origin
  outboxCoords$x<-outboxCoords$x*xyscale#convert pixels to real lengths
  outboxCoords$y<-outboxCoords$y*xyscale
  outboxCoords$z<-outboxCoords$z*zscale
  radius<-(((xydim*2)+1)/2)*xyscale
  distances<-as.data.frame(t(cdist(t(as.matrix(c(xydim*xyscale, xydim*xyscale, zdim*zscale))), as.matrix(outboxCoords))))
  if(nrow(distances)>0){#this prevents code from breaking if there are no pixels in vecinity
    distances<-subset(distances,!(distances[1]==0))#this removes the distance to self 
  }
  if(nrow(distances)>0){
    distances<-subset(distances,!(distances[1]>(radius)))#this filters out all the distances larger than the radius
  }
  #code crashes if you dont do independent "if statements"
  return(distances)
}

#The following 6 functionS are the third iteration of the pipeline. These functions have a spherical 
HollowSphereVolume<-function(r1,r2){
  volume1<-4/3*pi*(r1^3)
  volume2<-4/3*pi*(r2^3)
  return(volume2-volume1)
}
SpherePropOcc<-function(binsList, frequencyList, voxelSize){
  propOcc<-array(0,length(frequencyList))
  for (i in 1:length(frequencyList)){
    propOcc[i]<-frequencyList[i]*voxelSize/HollowSphereVolume(binsList[i],binsList[i+1])
  }
  return(propOcc)
}
RawBinningNorm<-function(dfDistances, boxVol, binSize, voxelVol){#takes the list of pairwise distances between focus and object in box and outputs frequency density table binned by distances of "binSize" in microns
  numDistances<-apply(dfDistances,1,as.numeric)#turn input df into numeric so hist() works
  bins<-seq(0, max(dfDistances)+binSize,binSize)#establish the binning criteria
  histogram<-hist(numDistances, breaks = bins, plot = FALSE)
  histCounts<-as.data.frame(histogram$counts)#these are the raw counts per bin
  histDensity<-as.data.frame(histogram$counts/length(numDistances))#these are the proportion of counts over non-empty pixels
  histDensityFunction<-as.data.frame(histogram$density*binSize)#these are the same as the top but using the function from hist class
  histProportion<-as.data.frame(histogram$counts/boxVol)#these are the proportion of counts over all pixels in box
  histNormSphereBinning<-as.data.frame(SpherePropOcc(bins,histogram$counts,voxelVol))
  finalData<-cbind(histogram$mids, histCounts, histDensity, histDensityFunction, histProportion, histNormSphereBinning)
  colnames(finalData)<- c('distance', 'counts', 'normDensity', "funcDensity", "rawDensity", "normSphereDensity")
  return(finalData)
}
NormRawLoopVoxels<-function(voxelsList, thresholdedImage, xyBoxDim, zBoxDim, xyScale, zScale, loopBinSize ){
  loopBox<-list()
  radius<-((xyBoxDim*2)+1)/2
  sphereVolum<-4/3*(radius^3)
  voxelVolum<-xyScale*xyScale*zScale
  for (i in 1:length(voxelsList[[1]])){
    distanceList<-DistanceRadius(voxelsList[i,],thresholdedImage, xyBoxDim, zBoxDim, xyScale, zScale)
    if(nrow(distanceList)>0){#this prevents code from breaking if there are no pixels in vecinity
      loopBox[[i]]<-RawBinningNorm(distanceList, sphereVolum,loopBinSize, voxelVolum)
    }
  }
  return(loopBox)
}
NormRawMergeStats<-function(dataList){
  RbindRawMergeLoopAll<-Reduce(function(dtf1,dtf2) rbind(dtf1,dtf2),dataList)
  RbindRawMergeLoopAll%>%group_by(distance)%>%summarise(meanCounts=mean(counts),normDensity=mean(normDensity),funcDensity=mean(funcDensity),rawDensity=mean(rawDensity), normPropOccup=mean(normSphereDensity))%>% as.data.frame()->finalData
  return(finalData)
}
NormRawPipeline<-function(redImg,greenImg,xysize,zsize,xydimSearch,zdimSearch,sampleSize,xyRealDim,zRealDim, pipeBinSize){
  voxCoordsRed<-CorrectSampleVolume(IdVoxel(redImg),xysize, zsize, xydimSearch, zdimSearch)
  voxCoordsGreen<-CorrectSampleVolume(IdVoxel(greenImg),xysize, zsize, xydimSearch, zdimSearch)
  print(paste0("Red pixels: ",nrow(voxCoordsRed)))
  print(paste0("Green pixels: ",nrow(voxCoordsGreen)))
  sampleVoxRed<-SampleVoxels(voxCoordsRed, sampleSize)
  sampleVoxGreen<-SampleVoxels(voxCoordsGreen, sampleSize)
  print("processing red to red")
  distListRR<-NormRawLoopVoxels(sampleVoxRed, redImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("processing red to green")
  distListRG<-NormRawLoopVoxels(sampleVoxRed, greenImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("processing green to red")
  distListGR<-NormRawLoopVoxels(sampleVoxGreen, redImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("processing green to green")
  distListGG<-NormRawLoopVoxels(sampleVoxGreen, greenImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  
  finalRR<-NormRawMergeStats(distListRR)
  finalRG<-NormRawMergeStats(distListRG)
  finalGR<-NormRawMergeStats(distListGR)
  finalGG<-NormRawMergeStats(distListGG)
  
  finalData<-gdata::combine(finalGG,finalRG,finalRR,finalGR)
  
  ggplot()+
    geom_point(data = finalData, aes(x=distance, y = rawDensity, color = source))+
    #geom_errorbar(data = finalData, aes(x=distance, ymin = mean-SD, ymax = mean+SD))
    ggtitle(paste("Density in search box r (XY=",((((xydimSearch*2)+1)/2)*xyRealDim)," Z=",((((zdimSearch*2)+1)/2)*zRealDim),")", sep = ''))
  
  ggplot()+
    geom_point(data = finalData, aes(x=distance, y = normPropOccup, color = source))+
    #geom_errorbar(data = finalData, aes(x=distance, ymin = mean-SD, ymax = mean+SD))
    ggtitle(paste("Density in search box r (XY=",((((xydimSearch*2)+1)/2)*xyRealDim)," Z=",((((zdimSearch*2)+1)/2)*zRealDim),")", sep = ''))
}
