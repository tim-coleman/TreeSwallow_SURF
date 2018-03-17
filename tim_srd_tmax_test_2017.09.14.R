
load("/Users/df36/projects/2017_Lucas_RF_CLT/data.surf/SRD_tmax_2014_200-365.RData")
str(dat)
load("/Users/df36/projects/2017_Lucas_RF_CLT/data.surf/surfing_eBird.srd.3km.locations.RData")
head(locations)

setwd("C:/Users/drain/OneDrive/eBird project")
load("C:/Users/drain/OneDrive/eBird project/ts.bcr30fall.RData")
set.seed(3000)
ts.perm <- ts.bcr30fall
ts.perm$temp.max <- sample(ts.bcr30fall$temp.max, replace = FALSE)

## Building the forests
library(randomForest)
rf_full  <- randomForest(occur~.-y-FID-lon-lat-date-dfs-year, data = ts.bcr30fall, ntree = 500)
rf_red <- randomForest(occur~.-y-FID-lon-lat-date-dfs-year, data = ts.perm, ntree = 500)

LID <- 1:length(dat[,1]) ## Creating location indices for easy subsetting. 
locations <- data.frame(locations, LID)

loc_NE <- locations[locations$lon > -78 & locations$lat > 37 & locations$lat < 44,]
srd_NE <- srd[[2]][loc_NE$LID,]
names(srd_NE) <- names(ts.bcr30fall)[11:26]



library(fields)
quilt.plot(x=loc_NE$lon, y=loc_NE$lat, z=dat2[,1])
map("state", add = T)
summary(dat[,50])

#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -30.0     4.5     7.5     7.9    11.5    27.5  687325

##########
# Noise
##########
FID <- rep(1, length(loc_NE$LID))
year <- FID
lon <- FID
lat <- FID
date <- FID
y <- FID
dfs <- FID


## Filling in the rest of the predictors

I.stationary <- rep(1, length(loc_NE$LID))
time <- rep(7.00, length(loc_NE$LID))
eff.hours <- rep(0.8405, length(loc_NE$LID))
eff.dist <- rep(1, length(loc_NE$LID))
no.obs <- rep(1, length(loc_NE$LID))

grid_df <- data.frame(srd_NE, I.stationary, time, eff.hours, eff.dist, no.obs,
                      FID, lon, lat, date, y, dfs, year)

## Putting in max temp information for various days
places <- loc_NE$LID
days <- seq(1,166,20)
dat2 <- dat[places, days]
names(dat2) <- paste(rep("Mtemp", 9), as.character(days), sep = "")


## Making the predictions
pred_diff <- function(tempmax){
  pred_full <- predict(rf_full, newdata = data.frame(grid_df, "temp.max" = tempmax), type = "response")
  pred_red <- predict(rf_red, newdata = data.frame(grid_df, "temp.max" = tempmax), type = "response")
  return(pred_full - pred_red)
}

## Making the predictions everywhere
pdiffs <- list()
for(i in 1:9){
  pdiff <- pred_diff(dat2[,i])
  pdiffs[[i]] <- data.frame("Pred_diff" = pdiff, "Lat" = loc_NE$lat, "Lon" = loc_NE$lon)
}

names(pdiffs) <- paste(rep("Mtemp", 9), as.character(days), sep = "")
save(pdiffs, file = "Mapping_Differences.RDA")

library(colorspace)
div_pal <- choose_palette()
dev.new()
levs <- sort(unique(c(seq(-.25,.25,.05), seq(-.09, .09, .01))))
par(mfrow = c(3,3), mar = c(4,3,3,4), oma = c(rep(.25,3), 1))
for(i in 1:9){
  quilt.plot(x=loc_NE$lon, y=loc_NE$lat, z=pdiffs[[i]]$Pred_diff,nlevel = length(levs)-1, breaks = levs,
             col = div_pal(length(levs)-1))
  title(paste("Day",days[i], "of Fall"))
  map("state", add = T)
}



