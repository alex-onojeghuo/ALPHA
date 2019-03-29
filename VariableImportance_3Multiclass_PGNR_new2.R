library(raster)
library(ggplot2)
library(rgdal)
library(dplyr)
library(plotly)
library(ggalt)
library(hrbrthemes)
library(caret)
library(randomForest)
library(SDMTools)
library(knitr)

rasterOptions(maxmemory = 1e+09,tmpdir = "X:/RtmpRasterDump")

lst <- c("VVdiff_summer_fall","NDWI_diff","B2_diff_summer_fall2","B3_diff_summer_fall2","B4_diff_summer_fall2","B8_diff_summer_fall2",
         "Topo_PC1","S2_PCA1", "S2_PCA2","S2_PCA3","S2_PCA4","Topo_PC4","TRI_LIDAR","SWI_LIDAR")
vars <- paste0(lst, ".tif")

OUTPUTS <- "X:/ALPHA/PGNR/RESULT"
setwd("X:/ALPHA/PGNR/Training")
p <- readOGR(".", "pts_new",  verbose=FALSE)
train <- raster("amwi_hr_traindata_final.tif")

##Location of input rasters for analysis
setwd("X:/5_PGNR_Mapping/1_Input_images_v2")
dat <- data.frame(row=1:length(p))
for (i in 1:length(vars)){
  r <- raster(vars[i])
  ext <- extract(r, p)
  dat <- cbind(dat,ext)
}
ext <- extract(train, p)
dat <- cbind(dat,ext)
dat <- dat[,-1]
d <- na.omit(dat)

joiner <- data.frame(WetlandTypeNum = c(1, 2, 3), WetlandType = c("Open Water", "Wetland", "Upland"), color = c('#006d2c', '#addd8e', '#fe9929'))
colnames(d) <- c(lst, "WetlandTypeNum")
df <- merge(d,joiner, by="WetlandTypeNum")

df.sub <- df[sample(1:length(df[,1]), size=2000, replace=FALSE),]


ctrl <-  trainControl(method = "repeatedcv",
                      repeats = 5,
                      number = 5)
RFModel <- train(as.factor(WetlandTypeNum)~., data = df.sub[,1:(length(lst)+1)], method = "rf", metric = "Accuracy", importance = T, ntree = 3000, trControl = ctrl)
ImpPlot <- plot(varImp(RFModel, scale = T)) #variable importance plots
ImpPlot
