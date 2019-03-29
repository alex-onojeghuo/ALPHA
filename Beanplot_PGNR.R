####################################################################
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#Bean plot chart preparation for data selection evaluation         #
#Filename: "InputDataExploration_Beanplot_V1.R"                    #
#Written and developed by Alex O. Onojeghuo                        #
#Alberta Biodiversity Monitoring Institute, Dec 1 2018             #
#------------------------------------------------------------------#
#------------------------------------------------------------------#
####################################################################

##Load libraries
library(ggplot2)
library(dplyr)
library(raster)
library(rgdal)
library(gtable)
library(grid)
library(zoo)
library(reshape2)
library(gridExtra)
library(beanplot)
##install.packages(c("raster", "sp"))

rasterOptions(maxmemory = 1e+09,tmpdir = "D:/RtmpRasterDump")

train.pts <- readOGR("X:/ALPHA/PGNR/Training","pts_new")

setwd("X:/5_PGNR_Mapping/1_Input_images_v2")

### List tif files in directory
list.files(pattern = "*.tif$")

###List required files
fls <- c("VVdiff_summer_fall.tif","NDWI_diff.tif","B2_diff_summer_fall2.tif","B3_diff_summer_fall2.tif","B4_diff_summer_fall2.tif",
         "B8_diff_summer_fall2.tif","Topo_PC1.tif","Topo_PC4.tif","TRI_LIDAR.tif","SWI_LIDAR.tif")

##Read in training data	
train <- raster("X:/ALPHA/PGNR/Training/amwi_hr_traindata_final.tif")

#Set output directory
OUTPUTS <- "X:/ALPHA/PGNR/RESULT"


##Extract training points from raster
df <- c(1:length(train.pts))
for (i in 1:length(fls)){
	r <- raster(fls[i])
	ext <- extract(r, train.pts)
	df <- cbind(df, ext)
	print(i)
}
##Extract training points from training data
ext <- extract(train, train.pts)
##Create dataframe between extracted points and training points
df <- cbind(df, ext)
##exclude first column from dataframe
df <- df[,-1]
##Remove all NA values
df <- na.omit(df)
##Define column names for dataframe
colnames(df) <- c("VV_diff","NDWI_diff","B2_diff","B3_diff","B4_diff","B8_diff","Topo_PC1","Topo_PC4","TRI","SWI","LC_class")

##If df is not a dataframe yet, convert to dataframe
df <- as.data.frame(df)
##check correlation of columns in dataframe
cor(df)
##Reclassify data with new class numbers


df$class <- ifelse(df$LC_class == 1, "Open Water", ifelse(df$LC_class ==2, "Wetland", "Upland"))

##Set working directory to store outputs
setwd(OUTPUTS)

##Split screen for plotting & output results as pdf
Funbean<-function(x){
X1<-cbind(df$class,x)
X1a<-X1[ which (X1$class == "Open Water"), ] 
X1b<-X1[ which (X1$class == "Wetland"), ]
X1c<-X1[ which (X1$class == "Upland"), ] 

X2<-as.data.frame(X1c[2,],X1b[,2],X1a[2,])}
X5<-list(df[,c(1,12)],df[c(2,12)],df[c(3,12)],df[c(4,12)],df[c(5,12)],df[c(6,12)],df[c(7,12)],df[c(8,12)],df[c(9,12)],df[c(10,12)])


names(X5)<-c("VV_diff","NDWI_diff","B2_diff","B3_diff","B4_diff","B8_diff","Topo_PC1","Topo_PC4","TRI","SWI")

#####################################
##1 VVdiff
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-10,10),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}

pdf(file="BeanPlot_New3.pdf",height=15 , width=10)
par(mfrow=c(3,4))

lapply(X5[1],Funbeana)

##2 NDWIdiff
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-0.5,0.5),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[2],Funbeana)

##3 B2_diff
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-1000,1000),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[3],Funbeana)

##4 B3_diff
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-1000,1000),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[4],Funbeana)

##5 B4_diff
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-1000,2000),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[5],Funbeana)

##6 B8_diff
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-3000,3000),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[6],Funbeana)


##7 Topo_PC1
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-5,5),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[7],Funbeana)

##8 Topo_PC4
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(-3,3),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[8],Funbeana)


##9 TRI
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(0,3),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}
lapply(X5[9],Funbeana)


##10 SWI
Funbeana<-function(X1){
  X1a<-X1[ which (X1$class == "Open Water"), ] 
  X1b<-X1[ which (X1$class == "Wetland"), ]
  X1c<-X1[ which (X1$class == "Upland"), ]
  beanplot(X1c[,1],X1b[,1], X1a[,1], main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" , ylim=c(2,15),
           log ="", col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50"), c("blanchedalmond", "gray50")), method="jitter")}


lapply(X5[10],Funbeana)


legend("topright",fill = c("darksalmon", "deepskyblue1"), legend = c("Open Water","Wetland","Upland"),bty = "n", cex=0.85)
dev.off()