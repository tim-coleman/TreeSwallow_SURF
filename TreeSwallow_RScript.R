####################################################################
# Welcome ! This is the R code associated with the paper:
# "Statistical Inferences on Tree Swallow Migrations, using Random Forests"
####################################################################

####################################################################
#         TABLE OF CONTENTS
####################################################################

#  0. Packages, helper functions
#  I. Introduction 
#     (1) Figure 1 Code
#  II. Data 
#     (1) Code associated with cleaning the Tree Swallow data.
#  III. Preliminary Models
#     (1) Figure 2 Code
#     (2) Figure 3 Code
#     (3) Figure 4 Code
#     (4) Moran's I testing Code
#  IV. Global Testing Procedure
#     (1) Figure 5 Code
#     (2) Figure 6 Code
#     (3) Global Testing Code
#  V. Local Testing Procedure
#     (1) Figure 7 Code
#     (2) Figure 8 Code
#     (3) Figure 9 Code
#     (4) Local Testing Code
#  VI. Discussion - No code here

####################################################################



## PACKAGES YOU WILL NEED ##
## Make sure you install each of these before loading

library(rpart)
library(randomForest)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

library(glmnet)
library(gam)
library(MASS)
library(FNN)
library(caret)
library(rpart)
library(Matrix)
library(brnn)
library(keras)
library(RSNNS)

library(rgdal)
library(sp)
library(rworldmap)
library(maps)
library(raster)

library(colorspace)
library(ape)
library(xtable)




# Helper Functions --------------------------------------------------------

## Color generating function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Helper function for subsetting
is.between <- function(x, a, b) {
  (x - a)  *  (b - x) > 0
}


####################################################################
## I INTRODUCTION
####################################################################


# (1) Figure 1 ------------------------------------------------------------
load('surfing_eBird.augmented.data.RData') #Loading the base dataset
# ----------------------------------------------------------------------
# BCR 30
# ----------------------------------------------------------------------	
png(file = paste(output.dir, "BCR.30.png",sep=""),
    height = 1000, 
    width = 1500)
erd.locs.index <- (erd$locs$bcr) == 30
sum(erd.locs.index, na.rm=T) 
# 435,038
map("world",xlim = c(-80, -65), 
    ylim = c(36,46), fill = TRUE, col = 'darkseagreen')
points(
  erd$locs$x[erd.locs.index], 
  erd$locs$y[erd.locs.index], 
  cex=0.5,col = 'yellow')
map("state", add=T)
title(main = "Bird Conservation Region 30")
# col="yellow",
# main = "Bird Conservation Region 30",
# xlab = "",
# ylab = "",
# add = T)
# -------
dev.off()



####################################################################
## II DATA
####################################################################
## This generates the basic BCR30 dataset.
ts.east <- data.frame("lon"=erd$locs$x[which(erd$locs$east.region==TRUE)],
                      "lat"=erd$locs$y[which(erd$locs$east.region==TRUE)],
                      "bcr"=erd$locs$bcr[which(erd$locs$east.region==TRUE)],
                      "date"=erd$dates[which(erd$locs$east.region==TRUE)],
                      "I.stationary"=erd$X$I.STATIONARY[which(erd$locs$east.region==TRUE)],
                      "year"=erd$X$YEAR[which(erd$locs$east.region==TRUE)],
                      "time"=erd$X$TIME[which(erd$locs$east.region==TRUE)],
                      "eff.hours"=erd$X$EFFORT_HRS[which(erd$locs$east.region==TRUE)],
                      "eff.dist"=erd$X$EFFORT_DISTANCE_KM[which(erd$locs$east.region==TRUE)],
                      "no.obs"=erd$X$NUMBER_OBSERVERS[which(erd$locs$east.region==TRUE)],
                      "aster.elev"=erd$X$ASTER2011_DEM[which(erd$locs$east.region==TRUE)],
                      "umd.LC"=erd$X$UMD2011_LANDCOVER[which(erd$locs$east.region==TRUE)],
                      "umd.C0"=erd$X$UMD2011_FS_C0_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C1"=erd$X$UMD2011_FS_C1_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C2"=erd$X$UMD2011_FS_C2_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C3"=erd$X$UMD2011_FS_C3_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C4"=erd$X$UMD2011_FS_C4_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C5"=erd$X$UMD2011_FS_C5_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C6"=erd$X$UMD2011_FS_C6_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C7"=erd$X$UMD2011_FS_C7_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C8"=erd$X$UMD2011_FS_C8_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C9"=erd$X$UMD2011_FS_C9_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C10"=erd$X$UMD2011_FS_C10_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C12"=erd$X$UMD2011_FS_C12_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C13"=erd$X$UMD2011_FS_C13_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "umd.C16"=erd$X$UMD2011_FS_C16_1500_PLAND[which(erd$locs$east.region==TRUE)],
                      "temp.min"=erd$X$tmin[which(erd$locs$east.region==TRUE)],
                      "temp.max"=erd$X$tmax[which(erd$locs$east.region==TRUE)],
                      "y"=erd$y[which(erd$locs$east.region==TRUE),3])


# Removing Missing Values ... 
# --------------------------------------------------------------
ts.east <- ts.east[complete.cases(ts.east),]

# --------------------------------------------------------------
# ... and Nonsensical Values
# --------------------------------------------------------------
ts.east <- ts.east[(ts.east$y>=0),]

# ------------------------------------------------------------------------
# Filtering to BCR 30, removing min temp and lat/lon
# ------------------------------------------------------------------------
ts.bcr30 <- ts.east[ts.east$bcr == 30,]
ts.bcr30 <- ts.bcr30[,-(which(names(ts.bcr30)=="temp.min" | names(ts.bcr30) == "bcr"))]

# ------------------------------------------------------------
# Filtering down to days after 200, getting occurrence
# ------------------------------------------------------------
frac.crit <- 200/365
date.frac <- ts.bcr30$date - ts.bcr30$year
indices <- which(date.frac >= frac.crit)
dfs <- date.frac[indices]
ts.bcr30fall <- ts.bcr30[indices,]
occur <- ifelse(ts.bcr30fall$y > 0, 1, 0) ## Creating occurence dataframe.
ts.bcr30fall <- cbind(ts.bcr30fall, occur)
ts.bcr30fall <- cbind(ts.bcr30fall, dfs)
FID <- 1:dim(ts.bcr30fall)[1]
ts.bcr30fall <- cbind(FID, ts.bcr30fall)
save(ts.bcr30fall, file = "ts.bcr30fall.RData")


####################################################################
## III PRELIMINARY MODELS
####################################################################


# Figure 2 Generation -----------------------------------------------------


library(randomForest)
library(glmnet)
library(gam)
library(MASS)
library(FNN)
library(caret)
library(dplyr)
library(rpart)
library(Matrix)
library(brnn)
library(keras)
library(RSNNS)

load("C:/Users/drain/OneDrive/eBird project/ts.bcr30fall.RData")

# Selecting only the predictors we want, and scaling them
ts.clean <- ts.bcr30fall %>% select(-c(FID, lat, lon, date, year, time, y))
ts.scale <- ts.clean %>% select(-occur) %>% mutate_all(funs(scale(.) %>% as.vector)) %>% 
  mutate(occur = ts.clean$occur)
write.csv(ts.scale, file = "birds.csv", row.names = F)

# Creating our 3 folds, and other training parameters.

train_control <- trainControl(method="cv", number=3, verboseIter = TRUE, returnData = FALSE)

# KNN
KNN_CV <- train(occur~., data = ts.scale, 
                method = "knn", tuneLength = 10, trControl = train_control)
save(KNN_CV, file = "KNN_CV.RDA")
min.KNN <- min(KNN_CV$results$RMSE)
# K = 11, RMSE = .2555827
load("C:/Users/drain/OneDrive/eBird project/KNN_CV.RDA")
# Random Forest
mtry <- sqrt(ncol(ts.scale)-1)
tg.rf <- expand.grid(.mtry=mtry)

RF_CV <- train(occur~., data = ts.scale,
               method = "rf", tuneGrid = tg.rf, trControl = train_control)
save(RF_CV, file = "RF_CV.RDA")

#mtry = 4.7, RMSE = .2227963

#GLM Net
tg.GLM = expand.grid(alpha = c(0,1),lambda = seq(0,0.15,by = 0.01))
GLMN_CV <- caret::train(occur~.*., data = ts.scale,
                        method = "glmnet", tuneGrid = tg.GLM, trControl = train_control)
save(GLMN_CV, file = "GLMN_CV.RDA")
# Lasso (alpha = 0), lambda = .01, RMSE =  0.2728308


#GAM
tg.GAM <- expand.grid(df = seq(1, 30, 5))
GAM_CV <- train(occur~., data = ts.scale,
                method = "gamSpline", tuneGrid = tg.GAM, trControl = train_control)
save(GAM_CV, file = "GAM_CV.RDA")
# df = 26, RMSE = .2738013

#Bayesian Reg. Neural Network
tg.BRNN <- expand.grid(neurons = c(5, 10, 30))
BRNN_CV <- train(occur~., data = ts.scale,
                 method = "brnn", tuneGrid = tg.BRNN, trControl = train_control)

#Multilayer Perceptron
tg.MP <- expand.grid(size = 2, rho = 100, dropout = .8, batch_size = 500, 
                     lr = .000001, decay = 0, activation = "relu")
MLP_CV <- train(as.factor(ifelse(occur==1, "B","NB"))~., data = ts.scale,
                method = "mlpKerasDropout", tuneGrid = tg.MP, trControl = train_class)
RMSE.calc(MLP_CV)


#Multilayer ANN
tg.ANN <- expand.grid(layer1 = c(15,30, 100), layer2 = c(15,30, 100), layer3 = c(15,30, 100))
ANN_CV <- caret::train(occur~., data = ts.scale,
                       method = "mlpML", tuneGrid = tg.ANN, trControl = train_control)
save(ANN_CV, file = "ANN_CV.RDA")

#RMSE = 0.2530037

## For methods that are classification only:

train_class <- trainControl(method="cv", number=3, verboseIter = TRUE, 
                            classProbs = TRUE, savePredictions = TRUE)
RMSE.calc <- function(mod){
  c("RMSE" = sqrt(mean((mod$pred$B - ts.scale$occur)^2)),
    "MAE" = mean(abs(mod$pred$B - ts.scale$occur)))
}

#LDA

LDA_CV <- train(as.factor(ifelse(occur==1, "B","NB"))~., data = ts.scale,
                method = "lda", trControl = train_class)
RMSE_LDA <- RMSE.calc(LDA_CV)
# RMSE_LDA = .3524737

# QDA
QDA_CV <- train(as.factor(ifelse(occur==1, "B","NB"))~., data = ts.scale,
                method = "qda", trControl = train_class)#QDA
RMSE_QDA <- RMSE.calc(QDA_CV)
# RMSE_QDA = .44427

library(xtable)
RF.r <- RF_CV$results[c("RMSE", "MAE")]
KNN.r <- KNN_CV$results[which.min(KNN_CV$results$RMSE),][c("RMSE", "MAE")]
ANN.r <- ANN_CV$results[which.min(ANN_CV$results$RMSE),][c("RMSE", "MAE")]
GAM.r <- GAM_CV$results[which.min(GAM_CV$results$RMSE),][c("RMSE", "MAE")]
GLM.r <- GLMN_CV$results[which.min(GLMN_CV$results$RMSE),][c("RMSE", "MAE")]

results.CV <- data.frame(rbind(RF.r, KNN.r, ANN.r,
                               GAM.r, GLM.r, RMSE_LDA,
                               RMSE_QDA))
row.names(results.CV) <- c("Random Forest \n mtry = 4.7", 
                           "KNN \n k = 11",
                           "ANN \n 100, 30, 15",
                           "GAM \n  df =26",
                           "L1 Logistic Reg., \n lambda = .01",
                           "LDA", "QDA")
models <- ordered(c("QDA", "LDA", "GAM", "GLMNet", "KNN", "ANN", "RF"))
models <- ordered(models, levels = c("QDA", "LDA", "GAM", "GLMNet", "KNN", "ANN", "RF"))

xtable(results.CV, digits = 5, align = "ccc")


library(ggplot2)
library(tibble)
library(tidyr)
results.CV.plot <- results.CV %>% arrange(desc(RMSE))%>% mutate(models = models) %>% 
  gather(Metric, Val, RMSE:MAE)

results.CV.plot %>% ggplot(aes(x = models, y = Val, fill = Metric)) +
  geom_col(position = "dodge", width = .5) + theme_classic() + ggtitle("3-Fold Cross Validation Error") +
  geom_hline(aes(yintercept = 0)) + coord_cartesian(ylim = c(0, .5)) + ylab("")+ xlab("")+
  theme(plot.title = element_text(size = 18, hjust = .5), axis.title = element_text(size = 15), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
        legend.background = element_rect(colour = "black", fill = "white"), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14), 
        legend.spacing.y = unit(2, "cm"), plot.subtitle = element_text(hjust = .5))


# Figure 3 ----------------------------------------------------------------


# -------------------------------------------------------------
# Generating the "ideal" test set
# -------------------------------------------------------------

date <- seq(frac.crit, 1, 1/365) #generates dates
tframe <- date
for(i in 2:23){ #assigns the mode of each predictor to our new test set
  index <- sample(1:2, 1)
  mval <- as.numeric(names(sort(table(ts.bcr30fall[,which(names(ts.bcr30fall) == names(ts.bcr30fall)[i])]), decreasing = TRUE)[index]))
  newvec <- rep(mval, length(date))
  print(newvec)
  tframe <- cbind(tframe, newvec)
}
tspline <-smooth.spline(ts.bcr30fall$temp.max~ I(ts.bcr30fall$date-ts.bcr30fall$year))
temps <- predict(tspline, date)
testset <- data.frame(tframe)
testset <- data.frame(tframe, temps$y)
names(testset) <- names(ts.bcr30fall)[1:24]
testset <- testset[,-which(names(testset) == "year")]

# Yearly Forests ----------------------------------------------------------

set.seed(1994)
train0809 <- trainset[trainset$year == 2008 | trainset$year == 2009,]
train1013 <- trainset[-which(trainset$year == 2008 | trainset$year == 2009),]
train1013_s <- train1013[sample(1:dim(train1013)[1], dim(train0809)[1]),]

trf.0809 <- randomForest(occur~. -year-lat-lon, data = train0809, ntree = 500)
trf.1013 <- randomForest(occur~. -year-lat-lon, data = train1013, ntree = 500)
trf.1013_s <- randomForest(occur~.-year-lat-lon, data = train1013_s, ntree = 500)

# -------------------------------------------------------------
# Generating the predictions on our test set
# -------------------------------------------------------------

pr0809 <- predict(rf.0809, newdata = testset, type = "response")
pr1013 <- predict(rf.1013, newdata = testset, type = "response")
pr1013_s <- predict(rf.1013_s, newdata = testset, type = "response")

plot(testset$date, pr1013, type = 's', ylim = c(0,150))
lines(testset$date, pr0809, col = 'red', type = 's')
lines(testset$date, pr1013_s, col = 'blue', type = 's')
text(.95, 50, "08-09 Curve", col = 'red')
text(.95, 45, "10-13 Curve", col = 'black')
text(.95, 40, "10-13 Curve, Reduced", col = 'blue')

pr.0809 <- predict(trf.0809, newdata = testset_2, type = "response")
pr.1013 <- predict(trf.1013, newdata = testset_2, type = "response")
pr.1013_s <- predict(trf.1013_s, newdata = testset_2, type = "response")

dfc <- cut(testset_2$dfs,date, labels = FALSE)
dzones <- cut(dfc, breaks = c(0,4,8,12), labels = FALSE)
tresult.0809 <- data.frame(dfc, pr.0809, dzones)[order(dfc),]
tresult.1013 <- data.frame(dfc, pr.1013, dzones)[order(dfc),]
tresult.1013_s <- data.frame(dfc, pr.1013_s, dzones)[order(dfc),]
t0809 <- tapply(tresult.0809$pr.0809, INDEX = dfc, FUN = mean)
t1013 <- tapply(tresult.1013$pr.1013, INDEX = dfc, FUN = mean)
t1013_s <- tapply(tresult.1013_s$pr.1013_s, INDEX = dfc, FUN = mean)
dts <- data.frame(t0809, t1013, t1013_s)
save(file = "RF_pred_pts", dts )


tresult.0809 <- data.frame(testset_2$dfs, pr.0809)[order(testset_2$dfs),]
tresult.1013 <- data.frame(testset_2$dfs, pr.1013)[order(testset_2$dfs),]
tresult.1013_s <- data.frame(testset_2$dfs, pr.1013_s)[order(testset_2$dfs),]
t0809 <- tapply(tresult.0809$pr.0809, INDEX = testset_2$dfs, FUN = mean)
t1013 <- tapply(tresult.1013$pr.1013, INDEX = testset_2$dfs, FUN = mean)
t1013_s <- tapply(tresult.1013_s$pr.1013_s, INDEX = testset_2$dfs, FUN = mean)
dts <- data.frame(t0809, t1013, t1013_s)


# Making the plots 

attach(dts)
smooth_0809 <- ksmooth(200:365, t0809, kernel = "normal", bandwidth = 10)$y
smooth_1013 <- ksmooth(200:365, t1013, kernel = "normal", bandwidth = 10)$y
smooth_1013_s <- ksmooth(200:365, t1013_s, kernel = "normal", bandwidth = 10)$y
detach(dts)


rf_labels <- as.factor(c(rep("2008-2009", 166), rep("2010-2013", 166), rep("2010-2013, reduced", 166)))
trf_df <- data.frame("Predictions" = c(smooth_0809, smooth_1013, smooth_1013_s), rf_labels)
rf_ps <- ggplot(data = trf_df, aes(x = rep(200:365,3), col = rf_labels, linetype = rf_labels)) + 
  geom_line(aes(y = Predictions), size = 1.25) +
  scale_x_continuous(name = "Day of Year", labels = seq(200, 365, 20), breaks = seq(200, 365,20)) +
  scale_y_continuous(name = "Occurrence", labels = seq(0, .3, .05), breaks = seq(0,.3,.05))+
  scale_color_manual(name = "Training Set", values =c(gg_color_hue(2), gg_color_hue(2)[2]),
                     labels = c("2008-2009", "2010-2013", "2010-2013\nreduced"))+
  scale_linetype_manual(name = "Training Set", values = c("solid", "dotted", "dashed"),
                        labels = c("2008-2009", "2010-2013", "2010-2013\nreduced")) + 
  geom_hline(aes(yintercept = 0)) + 
  ggtitle("Predicted Occurrence") + theme_classic()+
  theme(plot.title = element_text(size = 23, hjust = .5), axis.title = element_text(size = 17),
        axis.text = element_text(size = 13), legend.background = element_rect(colour = "black", fill = "white"),
        legend.title = element_text(size = 16), legend.spacing.y = unit(1.5, "cm"), 
        plot.subtitle = element_text(hjust = .5, size = 12), legend.text = element_text(size = 15), legend.key.size = unit(2.5,"line"))
rf_ps


# Figure 4 ----------------------------------------------------------------


ts.entry <- ts.bcr30fall[,-which(names(ts.bcr30fall) == "FID" | names(ts.bcr30fall) == "y"| 
                                   names(ts.bcr30fall) == "lon"  |names(ts.bcr30fall) == "lat"|
                                   names(ts.bcr30fall) == "date"|names(ts.bcr30fall) == "year")]
pr.occur <- partialPlot(rf.whole, ts.bcr30fall, x.var = temp.max)
plot(pr.occur, type = 'o', xlab = "Max Temp", 
     ylab = "Expected Occurence Conditional on Max Temp", main = "Partial Effect Plot for Max Temp")
abline(v = c(15, 17), lty = 2, col = 'red')
logits <- log(pr.occur$y/(1-pr.occur$y))
plot(pr.occur$x, logits, type = 'o', xlab = "Max Temp", ylab = "", main = "Partial Effect Plot for Max Temp")
save(rf.whole, file = "big_rf2.RDA")
save(pr.occur, file = "partial_effects.RDA")


load("partial_effects.RDA")
pr.occur <- data.frame(pr.occur)
ggplot(data = pr.occur, aes(x = pr.occur$x, y = pr.occur$y)) + geom_line(size = 2) + 
  xlab("Max Temp, Degrees Celsius") + ylab("Occurrence") + scale_x_continuous(breaks = seq(-15, 35,5)) + 
  ggtitle("Partial Effect of Maximum Temperature") +   theme_minimal()+
  theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
        axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15)) + 
  annotate(geom = "rect", xmin = 12.5, xmax = 17.5, ymin = .06, ymax = .1, alpha = .2, fill = 'red') +
  geom_curve(x = 0, xend = 12.5, y = .12, yend = .08,
             arrow = arrow(length = unit(0.3,"cm")), curvature = 0.5, size = .85)+
  geom_label(aes(x = 0, y = .12, label = 'Temperatures where insects\n begin to be active'), size = 5) 



# Figure 5 ----------------------------------------------------------------

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



# III (4) Moran's I Calculation -------------------------------------------


#Calculating Moran's I on each of the differenced plots
library(ape)

grid.dists <- as.matrix(dist(cbind(loc_NE$lon[seq(1, 34936, 6)], loc_NE$lat[seq(1, 34936, 6)])))

grid.dists.inv <- 1/grid.dists
diag(grid.dists.inv) <- 0

Istats <- list()
for(i in 1:9){
  Istats[[i]] <- Moran.I(pdiffs[[i]]$Pred_diff[seq(1, 34936, 6)], grid.dists.inv, na.rm = TRUE)
}

library(xtable)
I_mat <- unlist(Istats[[1]])
for(i in 2:9){
  I_row <- unlist(Istats[[i]])
  print(I_row)
  I_mat <-rbind(I_mat, I_row)
}
rownames(I_mat) <- NULL
I_mat <- data.frame(I_mat)
names(I_mat) <- c("Observed I", "Expected I", "Std Error", "P value")
Z_Score <- (I_mat$`Observed I`-I_mat$`Expected I`)/I_mat$`Std Error`
I_mat <- data.frame(I_mat, "Z" = Z_Score)
xtable(I_mat, digits = 4)

# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Observed.I & Expected.I & Std.Error & P.value & Z \\
# \hline
# 1 & 0.0422 & -0.0002 & 0.0004 & 0.0000 & 113.4010 \\
# 2 & 0.0407 & -0.0002 & 0.0004 & 0.0000 & 109.4391 \\
# 3 & 0.0302 & -0.0002 & 0.0004 & 0.0000 & 81.2810 \\
# 4 & 0.1045 & -0.0002 & 0.0004 & 0.0000 & 280.3697 \\
# 5 & 0.0757 & -0.0002 & 0.0004 & 0.0000 & 203.1577 \\
# 6 & 0.1673 & -0.0002 & 0.0004 & 0.0000 & 448.4178 \\
# 7 & 0.0724 & -0.0002 & 0.0004 & 0.0000 & 194.3888 \\
# 8 & 0.0638 & -0.0002 & 0.0004 & 0.0000 & 171.3538 \\
# 9 & 0.0804 & -0.0002 & 0.0004 & 0.0000 & 216.1006 \\
# \hline
# \end{tabular}
# \end{table}


####################################################################
## IV GLOBAL TESTING
####################################################################
#####
# Reducing the size of the 10-13 dataset, creating D^NN
####
ts.0809 <- ts.bcr30fall[which(ts.bcr30fall$year < 2010),]
ts.1013 <- ts.bcr30fall[which(ts.bcr30fall$year> 2009),]
lat08 <- scale(ts.0809$lat); lat10 <- scale(ts.1013$lat)
lon08 <- scale(ts.0809$lon); lon10 <- scale(ts.1013$lat)
dfs08 <- scale(ts.0809$dfs); dfs10 <- scale(ts.1013$dfs)
d08 <- data.frame("lat" = lat08, "lon" = lon08, "dfs" = dfs08)
d10 <- data.frame("lat" = lat10, "lon" = lon10, "dfs" = dfs10)

ts.NN_1013 <- ts.1013[50,]
for(i in 1:dim(ts.0809)[1]){
  obs0809 <- ts.0809[i,]
  ts.1013sub <- ts.1013[which(abs(ts.1013$dfs - obs0809$dfs) < .005 &
                                abs(ts.1013$lat - obs0809$lat) < .2 &
                                abs(ts.1013$lon - obs0809$lon) < .2),]
  if(dim(ts.1013sub)[1] < 1){
    ts.1013sub <- ts.1013[which(abs(ts.1013$dfs - obs0809$dfs) < .01 &
                                  abs(ts.1013$lat - obs0809$lat) < .25 &
                                  abs(ts.1013$lon - obs0809$lon) < .25),]
  }
  cat("\nNumber of observations triangulated: ", i)
  ts.1013NN <- ts.1013sub[sample(dim(ts.1013sub)[1], 1),]
  ts.NN_1013 <- rbind(ts.NN_1013, ts.1013NN)
}

save(ts.NN_1013, file = "C:/Users/drain/OneDrive/eBird project/ts.NN_1013.RDA")

#####
# Making partial effects by year and doing functional datatest
#####
load("C:/Users/drain/OneDrive/eBird project/ts.bcr30fall.RData")
ts.0809 <- ts.bcr30fall %>% filter(year < 2010)
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts.NN_1013.RDA")
load("C:/Users/drain/OneDrive/eBird project/SRD_tmax_1980-2007_200-365.RData")
load("C:/Users/drain/OneDrive/eBird project/surfing_eBird.srd.3km.locations.RData")
load("C:/Users/drain/OneDrive/eBird project/surfing_eBird.srd.3km.data.RData")
LID <- 1:length(dat[,1]) ## Creating location indices for easy subsetting. 
locations <- data.frame(locations, LID)
loc_NE <- locations[locations$lon > -78 & locations$lat > 37 & locations$lat < 44,]
srd_NE <- srd[[2]][loc_NE$LID,]
names(srd_NE) <- names(ts.bcr30fall)[11:26]
places <- loc_NE$LID
dat2 <- dat[places,]


# DoY model ---------------------------------------------------------------


test_gen_dfs <- function(nsample,day){
  samp_ind <- sample(1:34936, nsample, replace = FALSE)
  others <- srd_NE[samp_ind,]
  test <- data.frame(  "I.stationary" = rep(1, nsample*1),
                       "time" = rep(7.00, nsample*1),
                       "eff.hours" = rep(0.8405, nsample*1),
                       "eff.dist" = rep(1, nsample*1),
                       "no.obs" = rep(1, nsample*1),
                       others,
                       "dfs" = (day+199)/365)
  return(test[complete.cases(test),])
}


library(foreach)
library(randomForest)
library(dplyr)

set.seed(1994)
testpts <- list()
for(day in 1:166){
  testpts[[day]] <- test_gen_dfs(1000,day)
}

base_gen_dfs <- function(p, ntree, npar, test){
  junk <- c("FID", "year", "lon", "lat", "date", "y", "temp.max")
  ts.x0809 <- ts.0809[,-which(is.element(names(ts.0809), c(junk)))]
  ts.x1013 <- ts.NN_1013[,-which(is.element(names(ts.NN_1013), c(junk)))]
  # ts.y0809 <- ts.0809[,which(names(ts.0809) == "occur")]
  # ts.y0809 <- ts.NN_1013[,which(names(ts.NN_1013) == "occur")]
  ## Generating our test statistic
  start.time <- Sys.time()
  print(start.time)
  rf.1013 <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar% 
    randomForest(occur~., data = ts.x1013, ntree = ntree, mtry = p, nodesize = 15)
  rf.0809 <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar%  
    randomForest(occur~., data = ts.x0809, ntree = ntree, mtry = p, nodesize = 15)
  p0809.p <- lapply(testpts, predict, object = rf.1013)
  p1013.p <- lapply(testpts, predict, object = rf.0809)
  p0809.p <- unlist(lapply(p0809.p, FUN = mean))
  p1013.p <- unlist(lapply(p1013.p, FUN = mean))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(list("rf.1013" = rf.1013, "rf.0809" = rf.0809, "Preds0809" = p0809.p, "Preds1013" = p1013.p))
}

observed_RF_dfs <- base_gen_dfs(p = 22, ntree = 10, npar = 4)
save(observed_RF_dfs, file ="ObsRF_DFS.RDA")

pd_gen_dfs <- function(p, ntree, npar){
  junk <- c("FID", "year", "lon", "lat", "date", "y", "temp.max")
  start.time <- Sys.time()
  ts.bcr30fall_temp <- rbind(ts.0809, ts.NN_1013)
  ts.bcr30fall.p <- ts.bcr30fall_temp
  ts.bcr30fall.p$year <- sample(ts.bcr30fall_temp$year, replace = FALSE)
  ts.0809.p <- ts.bcr30fall.p[which(ts.bcr30fall.p$year < 2010),-which(is.element(names(ts.bcr30fall.p), junk))]
  ts.1013.p <- ts.bcr30fall.p[which(ts.bcr30fall.p$year> 2009),-which(is.element(names(ts.bcr30fall.p), junk))]
  # ts.y0809 <- ts.0809[,which(names(ts.0809) == "occur")]
  # ts.y0809 <- ts.NN_1013[,which(names(ts.NN_1013) == "occur")]
  ## Generating our test statistic
  start.time <- Sys.time()
  print(start.time)
  rf.1013 <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar% 
    randomForest(occur~., data = ts.1013.p, ntree = ntree, mtry = p, nodesize = 15)
  rf.0809 <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar%  
    randomForest(occur~., data = ts.0809.p, ntree = ntree, mtry = p, nodesize = 15)
  p0809.p <- lapply(testpts, predict, object = rf.1013)
  p1013.p <- lapply(testpts, predict, object = rf.0809)
  p0809.p <- unlist(lapply(p0809.p, FUN = mean))
  p1013.p <- unlist(lapply(p1013.p, FUN = mean))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(list("rf.1013_perm" = rf.1013, "rf.0809_perm" = rf.0809, 
              "Preds0809_perm" = p0809.p, "Preds1013_perm" = p1013.p))
}

perm_RF_dfs <- pd_gen_dfs(p = 22, ntree = 10, npar = 4)



set.seed(19943)
for(i in 1:1000){
  setwd("~/PermutatedDFS")
  pd1 <- pd_gen_dfs(p = 22, ntree = 10, npar = 4)
  fname <- paste("~/PermutatedDFS/Permutation_", i, "_DFS.RDA",sep = "")
  save(pd1, file = fname)
}

## Creating the prediction dataframes.
setwd("~/PermutatedDFS")
dfs_dir <- dir()
p0809_dfs <- rep(0, 166)
p1013_dfs <- rep(0, 166)
for(i in 1:length(dfs_dir)){
  load(dfs_dir[i])
  p0809s <- lapply(testpts, predict, object = pd1$rf.0809)
  p1013s <- lapply(testpts, predict, object = pd1$rf.1013)
  p0809s <- unlist(lapply(p0809s, FUN = mean))
  p1013s <- unlist(lapply(p1013s, FUN = mean))
  p0809_dfs <- cbind(p0809_dfs, p0809s)
  p1013_dfs <- cbind(p1013_dfs, p1013s)
}

p0809.dfs <- p0809_dfs[,-1]
p1013.dfs <- p1013_dfs[,-1]

setwd("C:/Users/drain/OneDrive/eBird project")
save(object = p0809.dfs, file = "preds_0809dfs.RData")
save(object = p1013.dfs, file ="preds_1013dfs.RData")



# Full (DoY & Max Temp) Model Testing ------------------------------------------------------



test_gen_full <-function(nsample,day){
  samp_ind <- sample(1:34936, nsample, replace = FALSE)
  others <- srd_NE[samp_ind,]
  temps <- dat2[samp_ind, day]
  others <- srd_NE[samp_ind,]
  test <- data.frame(  "I.stationary" = rep(1, nsample*1),
                       "time" = rep(7.00, nsample*1),
                       "eff.hours" = rep(0.8405, nsample*1),
                       "eff.dist" = rep(1, nsample*1),
                       "no.obs" = rep(1, nsample*1),
                       others,
                       "temp.max" = temps,
                       "dfs" = (day+199)/365)
  return(test[complete.cases(test),])
}

set.seed(1994)
testpts <- list()
for(day in 1:166){
  testpts[[day]] <- test_gen_full(1000,day)
}
setwd("C:/Users/drain/OneDrive/eBird project")
observed_RF_full <- base_gen_dfs(p = 23, ntree = 10, npar = 4)
save(observed_RF_full, file ="ObsRF_full.RDA")


for(i in 1:1000){
  setwd("~/Full_Permuted")
  pd1 <- pd_gen_dfs(p = 23, ntree = 10, npar = 4)
  fname <- paste("~/Full_Permuted/Permutation_", i, ".RDA",sep = "")
  save(pd1, file = fname)
}


## Creating the prediction dataframes.
setwd("~/Full_Permuted")
dfs_dir <- dir()
p0809_dfs <- rep(0, 166)
p1013_dfs <- rep(0, 166)
for(i in 1:length(dfs_dir)){
  load(dfs_dir[i])
  p0809s <- lapply(testpts, predict, object = pd1$rf.0809)
  p1013s <- lapply(testpts, predict, object = pd1$rf.1013)
  p0809s <- unlist(lapply(p0809s, FUN = mean))
  p1013s <- unlist(lapply(p1013s, FUN = mean))
  p0809_dfs <- cbind(p0809_dfs, p0809s)
  p1013_dfs <- cbind(p1013_dfs, p1013s)
  if(i %in% seq(0,1000,50)){print(i)}
}

p0809.full <- p0809_dfs[,-1]
p1013.full <- p1013_dfs[,-1]

setwd("C:/Users/drain/OneDrive/eBird project")
save(object = p0809.full, file = "preds_0809_full.RData")
save(object = p1013.full, file ="preds_1013_full.RData")



pd_gen <- function(nperm, p, ntree, npar){
  
  ## Generating our idealized testing set
  places <- loc_NE$LID
  dat2 <- dat[places,]
  test.ideal <- data.frame(
    "temp.max" = apply(dat2, FUN = mean, MARGIN = 2, na.rm = TRUE),
    "I.stationary" = rep(1, 166),
    "time" = rep(7.00, 166),
    "eff.hours" = rep(0.8405, 166),
    "eff.dist" = rep(1, 166),
    "no.obs" = rep(1, 166),
    "aster.elev" = mean(srd_NE$aster.elev),
    "umd.LC"= mean(srd_NE$umd.LC),
    "umd.C0" = mean(srd_NE$umd.C0),
    "umd.C1" = mean(srd_NE$umd.C1),
    "umd.C2" = mean(srd_NE$umd.C2),
    "umd.C3" = mean(srd_NE$umd.C3),
    "umd.C4" = mean(srd_NE$umd.C4),
    "umd.C5" = mean(srd_NE$umd.C5),
    "umd.C6" = mean(srd_NE$umd.C6),
    "umd.C7" = mean(srd_NE$umd.C7),
    "umd.C8" = mean(srd_NE$umd.C8),
    "umd.C9" = mean(srd_NE$umd.C9),
    "umd.C10" = mean(srd_NE$umd.C10),
    "umd.C12" = mean(srd_NE$umd.C12),
    "umd.C13" = mean(srd_NE$umd.C13),
    "umd.C16" = mean(srd_NE$umd.C16))
  junk <- c("FID", "year", "lon", "lat", "date", "y", "dfs")
  ## Permutation Distribution
  CvMs <- c(); KSs <- c()
  perm0809 <- matrix(nrow = 166, ncol = nperm)
  perm1013 <- matrix(nrow = 166, ncol = nperm)
  for(perm in 1:nperm){
    start.time <- Sys.time()
    ts.bcr30fall_temp <- rbind(ts.0809, ts.NN_1013)
    ts.bcr30fall.p <- ts.bcr30fall_temp
    ts.bcr30fall.p$year <- sample(ts.bcr30fall_temp$year, replace = FALSE)
    ts.0809.p <- ts.bcr30fall.p[which(ts.bcr30fall.p$year < 2010),-which(is.element(names(ts.bcr30fall.p), junk))]
    ts.1013.p <- ts.bcr30fall.p[which(ts.bcr30fall.p$year> 2009),-which(is.element(names(ts.bcr30fall.p), junk))]
    ts.x0809.p <- ts.0809.p[,-which(names(ts.0809.p) == "occur")]
    ts.x1013.p <- ts.1013.p[,-which(names(ts.1013.p) == "occur")]
    ts.y0809.p <- ts.0809.p[,which(names(ts.0809.p) == "occur")]
    ts.y1013.p <- ts.1013.p[,which(names(ts.1013.p) == "occur")]
    ## Generating our test statistic
    rf.1013.p <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar% 
      randomForest(x = ts.x1013.p, y = ts.y1013.p,ntree = ntree, mtry = p, nodesize = 15)
    rf.0809.p <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar%  
      randomForest(x = ts.x0809.p, y = ts.y0809.p,ntree = ntree, mtry = p, nodesize = 15)
    pred.0809.p <- predict(rf.0809.p, newdata = test.ideal)
    pred.1013.p <- predict(rf.1013.p, newdata = test.ideal)
    KSs[perm] <- max(abs(pred.0809.p-pred.1013.p))
    CvMs[perm] <- mean((pred.0809.p-pred.1013.p)^2)
    perm0809[,nperm] <- pred.0809.p
    perm1013[,nperm] <- pred.1013.p
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }
  return(list("Permuted CvMs" = CvMs, "Permuted KS" = KSs,
              "Permuted0809" = perm0809, "Permuted1013" = perm1013, 
              "RF.0809" <- rf.0809.p, "RF.1013" <- rf.1013.p))
}

base_gen_mtemp <- function(p, ntree, npar){
  places <- loc_NE$LID
  dat2 <- dat[places,]
  test.ideal <- data.frame(
    "temp.max" = apply(dat2, FUN = mean, MARGIN = 2, na.rm = TRUE),
    "I.stationary" = rep(1, 166),
    "time" = rep(7.00, 166),
    "eff.hours" = rep(0.8405, 166),
    "eff.dist" = rep(1, 166),
    "no.obs" = rep(1, 166),
    "aster.elev" = mean(srd_NE$aster.elev),
    "umd.LC"= mean(srd_NE$umd.LC),
    "umd.C0" = mean(srd_NE$umd.C0),
    "umd.C1" = mean(srd_NE$umd.C1),
    "umd.C2" = mean(srd_NE$umd.C2),
    "umd.C3" = mean(srd_NE$umd.C3),
    "umd.C4" = mean(srd_NE$umd.C4),
    "umd.C5" = mean(srd_NE$umd.C5),
    "umd.C6" = mean(srd_NE$umd.C6),
    "umd.C7" = mean(srd_NE$umd.C7),
    "umd.C8" = mean(srd_NE$umd.C8),
    "umd.C9" = mean(srd_NE$umd.C9),
    "umd.C10" = mean(srd_NE$umd.C10),
    "umd.C12" = mean(srd_NE$umd.C12),
    "umd.C13" = mean(srd_NE$umd.C13),
    "umd.C16" = mean(srd_NE$umd.C16)
    # "FID"= rep(1, 166),
    # "year" = rep(1,166),
    # "lon" = rep(1,166), ## Junk in the model
    # "lat" = rep(1,166),
    # "date" = rep(1,166),
    # "y" = rep(1,166),
    # "dfs" = rep(1,166)
  )
  junk <- c("FID", "year", "lon", "lat", "date", "y", "dfs")
  ts.bcr30fall_temp <- rbind(ts.0809, ts.NN_1013)
  ts.0809 <- ts.bcr30fall_temp[which(ts.bcr30fall_temp$year < 2010),-which(is.element(names(ts.bcr30fall_temp), junk))]
  ts.1013 <- ts.bcr30fall[which(ts.bcr30fall_temp$year> 2009),-which(is.element(names(ts.bcr30fall_temp), junk))]
  ts.x0809 <- ts.0809[,-which(names(ts.0809) == "occur")]
  ts.x1013 <- ts.1013[,-which(names(ts.1013) == "occur")]
  ts.y0809 <- ts.0809[,which(names(ts.0809) == "occur")]
  ts.y1013 <- ts.1013[,which(names(ts.1013) == "occur")]
  ## Generating our test statistic
  start.time <- Sys.time()
  print(start.time)
  rf.1013 <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar% 
    randomForest(x = ts.x1013, y = ts.y1013,ntree = ntree, mtry = p, nodesize = 15)
  rf.0809 <- foreach(ntree=rep(ntree, npar), .combine=randomForest::combine, .packages='randomForest') %dopar%  
    randomForest(x = ts.x0809, y = ts.y0809,ntree = ntree, mtry = p, nodesize = 15)
  pred.0809 <- predict(rf.0809, newdata = test.ideal)
  pred.1013 <- predict(rf.1013, newdata = test.ideal)
  KS.obs <- max(abs(pred.0809-pred.1013))
  CvM.obs <- mean((pred.0809-pred.1013)^2)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(list("rf1013" = rf.1013, "rf.0809" = rf.0809, "Preds0809" = pred.0809, "Preds1013" = pred.1013, "Observed KS Statistic" = KS.obs,
              "Observed CvM statistic" = CvM.obs))
}

observed_RF <- base_gen_mtemp(p = 22, ntree = 10, npar = 4)
save(object = observed_RF, file = "observed_RF.Rdata")



set.seed(19943)
perm_ind <- list()
for(i in 439:1000){
  setwd("C:/Users/drain/OneDrive/eBird project/Permutations")
  pd1 <- pd_gen(nperm = 1, p = 22, ntree = 10, npar = 4)
  fname <- paste("C:/Users/drain/OneDrive/eBird project/Permutations/Permutation_", i, ".RDA",sep = "")
  save(pd1, file = fname)
}



permtest_1000 <- pd_gen(1000)
setwd("C:/Users/drain/OneDrive/eBird project/Permutations")
d1 <- dir()
KS_p <- c(); CvM_p <- c()
for(i in 1:length(d1)){
  load(d1[i])
  KS_p[i] <- pd1$`Permuted KS`; CvM_p[i] <- pd1$`Permuted CvMs`
}


## not so good ^^, so we generate another idealized test set that is more representative of the observations.
load("C:/Users/drain/OneDrive/eBird project/SRD_tmax_1980-2007_200-365.RData")
load("C:/Users/drain/OneDrive/eBird project/surfing_eBird.srd.3km.locations.RData")
load("C:/Users/drain/OneDrive/eBird project/surfing_eBird.srd.3km.data.RData")
load("C:/Users/drain/OneDrive/eBird project/observed_RF.Rdata")

LID <- 1:length(dat[,1]) ## Creating location indices for easy subsetting. 
locations <- data.frame(locations, LID)
loc_NE <- locations[locations$lon > -78 & locations$lat > 37 & locations$lat < 44,]
places <- loc_NE$LID
srd_NE <- srd[[2]][loc_NE$LID,]
names(srd_NE) <- names(ts.bcr30fall)[11:26]


dat2 <- dat[places, 1:166]

test_gen <- function(nsample,day){
  samp_ind <- sample(1:34936, nsample, replace = FALSE)
  temps <- dat2[samp_ind, day]
  others <- srd_NE[samp_ind,]
  test <- data.frame(  "I.stationary" = rep(1, nsample*1),
                       "time" = rep(7.00, nsample*1),
                       "eff.hours" = rep(0.8405, nsample*1),
                       "eff.dist" = rep(1, nsample*1),
                       "no.obs" = rep(1, nsample*1),
                       others,
                       "temp.max" = temps)
  return(test[complete.cases(test),])
}

set.seed(1994)
testpts <- list()
for(day in 1:166){
  testpts[[day]] <- test_gen(1000,day)
}
library(randomForest)
p0809 <- lapply(testpts, predict, object = observed_RF$rf.0809)
p1013 <- lapply(testpts, predict, object = observed_RF$rf1013)
p0809 <- unlist(lapply(p0809, FUN = mean))
p1013 <- unlist(lapply(p1013, FUN = mean))
save(p0809, file = "p0809.RDA")
save(p1013, file = "p1013.RDA")
ks_obs <- max(abs(p0809 - p1013))
cvm_obs <- mean((p0809-p1013)^2)
r_pobs <- cor(p0809, p1013, method = "pearson")
r_sobs <- cor(p0809, p1013, method = "spearman")

plot.ts(ksmooth(x = 1:166, p0809, kernel = "normal", bandwidth = 5)$y)
lines(ksmooth(x = 1:166, p1013, kernel = "normal", bandwidth = 5)$y, col ='red')

df_preds <- data.frame("preds" = 
                         c(ksmooth(x = 1:166, p0809, kernel = "normal", bandwidth =7.5)$y, 
                           ksmooth(x = 1:166, p1013, kernel = "normal", bandwidth = 7.5)$y), "yrs" <- c(rep("0809", 166), rep("1013", 166)))
gpreds <- ggplot(data = df_preds, aes(x = rep(1:166,2), y = preds, color = yrs)) + geom_line(lwd = 1.15) +
  scale_color_discrete(name = "Training\nYears", labels = c("2008-2009","2010-2013")) + 
  scale_x_continuous(name = "Day of Year", breaks = seq(0, 165, 15)) + 
  scale_y_continuous(name = "Occurrence", breaks = seq(0, .3, .05)) + ggtitle("Observed Prediction Curves") + theme_minimal()+
  theme(plot.title = element_text(size = 16, hjust = .5), axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), legend.background = element_rect(colour = "black", fill = "white"),
        legend.title = element_text(size = 13), legend.spacing.y = unit(1, "cm"), legend.text = element_text(size = 12)) 
gpreds

d1 <- dir()
ks_p <- c(); cvm_p <- c(); r_p <- c(); r_s <- c()
preds_0809 <- rep(0, 166)
preds_1013 <- rep(0, 166)
for(i in 1:length(d1)){
  load(d1[i])
  p0809.p <- lapply(testpts, predict, object = pd1[[5]])
  p1013.p <- lapply(testpts, predict, object = pd1[[6]])
  p0809.p <- unlist(lapply(p0809.p, FUN = mean))
  p1013.p <- unlist(lapply(p1013.p, FUN = mean))
  preds_0809 <- cbind(preds_0809, p0809.p) 
  preds_1013 <- cbind(preds_1013, p1013.p)
  ks_p[i] <- max(abs(p0809.p - p1013.p))
  cvm_p[i] <- mean((p0809.p-p1013.p)^2)
  r_p[i] <- cor(p0809.p, p1013.p, method = "pearson")
  r_s[i] <- cor(p0809.p, p1013.p, method = "spearman")
}

preds_0809 <- preds_0809[,-1]
preds_1013 <- preds_1013[,-1]
save(object = preds_0809, file = "preds_0809.RData")
save(object = preds_1013, file ="preds_1013.RData")

## smoothed_stats is a function that takes in dataframes, returns p-values, plots of estimated functions, and histograms.

smoothed_stats <- function(dfo1, dfo2, dfp1, dfp2, bw, index, max.y = .45, unsmooth = FALSE, smooth.stats = F){
  if(unsmooth){
    cvm_smooth <- c(); ks_smooth <- c(); D1 <- c(); D2 <- c(); D3 <- c()
    for(i in 1:dim(dfp1)[2]){
      D1[i] <- mean(dfp1[1:55,i] - dfp2[1:55,i])
      D2[i] <- mean(dfp1[56:110,i] - dfp2[56:110,i])
      D3[i] <- mean(dfp1[111:166,i] - dfp2[111:166,i])
      ks_smooth[i] <- max(abs(dfp1[,i] - dfp2[,i]))
      cvm_smooth[i] <- mean((dfp1[,i] - dfp2[,i])^2)
    }
    cvm_obs = mean((dfo2 - dfo1)^2); ks_obs = max(abs(dfo2 - dfo1))
    D1_obs <- mean(dfo1[1:65] - dfo2[1:65])
    D2_obs <- mean(dfo1[66:110] - dfo2[66:110])
    D3_obs <- mean(dfo1[111:166] - dfo2[111:166])
    df_tstat <- data.frame(CvM =  50*cvm_smooth, KS = ks_smooth)
    df_tstat <- melt(df_tstat)
    gdensityCvMKS <- ggplot(data = df_tstat, aes(x = value, fill = variable)) + 
      geom_density(mapping = aes(linetype = variable), alpha = .45, size = 1.095) +
      scale_fill_manual(values=c("blue", "Tomato"),name = "Test Statistic", 
                        labels = c("Cramer von Mises\n", "Kolmorgorov-\nSmirnov")) +
      scale_linetype_manual(values=c("solid", "dotdash"),name = "Test Statistic", 
                            labels = c("Cramer von Mises\n", "Kolmorgorov-\nSmirnov")) +
      geom_segment(aes(x = ks_obs, y = 0, xend =ks_obs, yend = 8), lwd = 1.85, col = "Tomato", linetype = "dotdash")+ 
      geom_segment(aes(x = 50*cvm_obs, y = 0, xend =50*cvm_obs, yend = 3), lwd = 1.85, col = "blue", linetype = "solid") + 
      coord_cartesian(xlim = c(0, .15), ylim = c(0, 25)) + ylab(NULL)+
      theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
            axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 15),
            axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15),
            legend.background = element_rect(colour = "black", fill = "white"),
            legend.title = element_text(size = 15), legend.spacing.y = unit(.75, "cm"),
            legend.text = element_text(size = 13), legend.key.size = unit(2,"line"))  +
      scale_x_continuous(name = "Test Statistic") + ggtitle("Estimated Permutation Densities")
    df_dstat <- melt(data.frame(D1, D2, D3))
    gD2<- ggplot(data = df_dstat, aes(x = value, fill = variable), size = 1.5) +
      geom_density(mapping = aes(linetype = variable), alpha = .45, size = 1.095) +
      scale_fill_manual(values=rev(c("#002984", "#D240A7" ,"#FFE072")),name = "Test Statistic", 
                        labels = c("Curve Distance,\n DoY 200-264\n", 
                                   "Curve Distance,\n DoY 265-310\n", "Curve Distance,\n DoY 310-365\n")) +
      scale_linetype_manual(values=c("twodash", "solid", "dashed"),name = "Test Statistic", 
                            labels = c("Curve Distance,\n DoY 200-264\n", 
                                       "Curve Distance,\n DoY 265-310\n", "Curve Distance,\n DoY 310-365\n")) +
      geom_segment(aes(x = D1_obs, y = 0, xend =D1_obs, yend = 6), lwd = 1.85, col = "#FFE072", linetype = "twodash")+ 
      geom_segment(aes(x = D2_obs, y = 0, xend =D2_obs, yend = 2.5), lwd = 1.85, col = "#D240A7", linetype = "solid") + 
      geom_segment(aes(x = D3_obs, y = 0, xend =D3_obs, yend = 19), lwd = 1.85, col = "#002984", linetype = "dashed") + 
      coord_cartesian(xlim = c(-.15, .15), ylim = c(0, 25)) + ylab(NULL)+
      theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
            axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 15),
            axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15),
            legend.background = element_rect(colour = "black", fill = "white"),
            legend.title = element_text(size = 15), legend.spacing.y = unit(.75, "cm"),
            legend.text = element_text(size = 13), legend.key.size = unit(2,"line"))  +
      scale_x_continuous(name = "Test Statistic") + ggtitle("Estimated Permutation Densities")
    df_preds <- data.frame("preds" = c(dfo1, dfo2), "yrs" <- c(rep("0809", 166), rep("1013", 166)))
    gpreds <- ggplot(data = df_preds, aes(x = rep(200:365,2), y = preds, color = yrs)) + geom_line(lwd = 1.15) +
      scale_color_discrete(name = "Training\nYears", labels = c("2008-2009","2010-2013")) + 
      coord_cartesian(xlim = c(200, 365), ylim = c(0, max.y)) +
      scale_x_continuous(name = "Day of Year", breaks = seq(200, 365, 15)) + 
      scale_y_continuous(name = "Occurrence", breaks = seq(0, max.y, .05)) + ggtitle("Observed Prediction Curves") + 
      theme_classic()+
      theme(plot.title = element_text(size = 20, hjust = .5), axis.title = element_text(size = 17),
            axis.text = element_text(size = 17), legend.background = element_rect(colour = "black", fill = "white"),
            legend.title = element_text(size = 17), legend.spacing.y = unit(1, "cm"),
            legend.text = element_text(size = 15),legend.key.size = unit(2,"line")) 
    df_perm <- data.frame("preds" = c(dfp1[,index], dfp2[,index]),"yrs" = c(rep("0809", 166), rep("1013", 166)))
    gperm <- ggplot(data = df_perm, aes(x = rep(200:365,2), y = preds, color = yrs)) + geom_line(lwd = 1.15) +
      scale_color_discrete(name = "Training\nYears", labels = c("2008-2009","2010-2013")) + 
      coord_cartesian(xlim = c(200, 365), ylim = c(0, max.y)) +
      scale_x_continuous(name = "Day of Year", breaks = seq(200, 365, 15)) + 
      scale_y_continuous(name = "Occurrence", breaks = seq(0, max.y, .05)) + ggtitle("Permuted Prediction Curves") +
      theme_classic()+
      theme(plot.title = element_text(size = 20, hjust = .5), axis.title = element_text(size = 17),
            axis.text = element_text(size = 17), legend.background = element_rect(colour = "black", fill = "white"),
            legend.title = element_text(size = 17), legend.spacing.y = unit(1, "cm"),
            legend.text = element_text(size = 15),legend.key.size = unit(2,"line")) 
    return(list(perm_dist = data.frame("CvM" = cvm_smooth, "KS" = ks_smooth, "D1" = D1, "D2" = D2, "D3" = D3),
                "obs_stats" = c("CvM" = cvm_obs, "KS" = ks_obs, "Diff_1" = D1_obs, "Diff_2" = D2_obs, "Diff_3" = D3_obs),
                "Density_Plot" = gdensityCvMKS,  "Predicted_Plot" = gpreds, "Permuted_Plot" = gperm, "D1to3Density" = gD2))
  }
  else{
    dfo1 <- ksmooth(1:166, dfo1, kernel = 'normal', bandwidth = bw)$y
    dfo2 <- ksmooth(1:166,dfo2, kernel = 'normal', bandwidth = bw)$y
    dfp1 <- apply(dfp1, FUN = ksmooth, MARGIN = 2, kernel = 'normal', bandwidth = bw, x=1:166)
    dfp2 <- apply(dfp2, FUN = ksmooth, MARGIN = 2, kernel = 'normal', bandwidth = bw, x=1:166)
    if (smooth.stats) {
      cvm_smooth <- c(); ks_smooth <- c(); D1 <- c(); D2 <- c(); D3 <- c()
      for(i in 1:length(dfp1)){
        ks_smooth[i] <- max(abs(dfp1[i]$p0809.p$y - dfp2[i]$p1013.p$y))
        cvm_smooth[i] <- mean((dfp1[i]$p0809.p$y - dfp2[i]$p1013.p$y)^2)
        D1[i] <- mean(dfp1[i]$p0809.p$y[1:55] - dfp2[i]$p1013.p$y[1:55])
        D2[i] <- mean(dfp1[i]$p0809.p$y[56:110] - dfp2[i]$p1013.p$y[56:110])
        D3[i] <- mean(dfp1[i]$p0809.p$y[111:166] - dfp2[i]$p1013.p$y[111:166])
      }
      cvm_obs = mean((dfo2 - dfo1)^2); ks_obs = max(abs(dfo2 - dfo1))
      D1_obs <- mean(dfo1[1:65] - dfo2[1:65])
      D2_obs <- mean(dfo1[66:110] - dfo2[66:110])
      D3_obs <- mean(dfo1[111:166] - dfo2[111:166])
      df_tstat <- data.frame(CvM =  50*cvm_smooth, KS = ks_smooth)
      df_tstat <- melt(df_tstat)
      gdensity <- ggplot(data = df_tstat, aes(x = value, fill = variable), lwd = 0) + geom_density(alpha = .5) +
        scale_fill_manual(values=c("blue", "Tomato"),name = "Test Statistic", 
                          labels = c("Cramer von Mises\n", "Kolmorgorov-\nSmirnov")) +
        geom_segment(aes(x = ks_obs, y = 0, xend =ks_obs, yend = 5), lwd = 1.5, col = "Tomato")+ 
        geom_segment(aes(x = 50*cvm_obs, y = 0, xend =50*cvm_obs, yend = 30), lwd = 1.5, col = "blue") + 
        coord_cartesian(xlim = c(0, .15), ylim = c(0, 30)) + ylab(NULL)+
        theme_minimal() +
        theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
              axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 15),
              axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15),
              legend.background = element_rect(colour = "black", fill = "white"),
              legend.title = element_text(size = 15), legend.spacing.y = unit(.75, "cm"),
              legend.text = element_text(size = 13), legend.key.size = unit(1.5,"line"))  +
        scale_x_continuous(name = "Test Statistic") + ggtitle("Estimated Permutation Densities")
      df_dstat <- melt(data.frame(D1, D2, D3))
      gD2<- ggplot(data = df_dstat, aes(x = value, fill = variable), lwd = 0) + geom_density(alpha = .5) +
        scale_fill_manual(values=rev(c("#002984", "#D240A7" ,"#FFE072")),name = "Test Statistic", 
                          labels = c("Curve Distance,\n DoY 200-264\n", 
                                     "Curve Distance,\n DoY 265-310\n", "Curve Distance,\n DoY 310-365\n")) +
        geom_segment(aes(x = D1_obs, y = 0, xend =D1_obs, yend = 45), lwd = 1.5, col = "#FFE072")+ 
        geom_segment(aes(x = D2_obs, y = 0, xend =D2_obs, yend = 45), lwd = 1.5, col = "#D240A7") + 
        geom_segment(aes(x = D3_obs, y = 0, xend =D3_obs, yend = 45), lwd = 1.5, col = "#002984") + 
        coord_cartesian(xlim = c(-.15, .15), ylim = c(0, 25)) + ylab(NULL)+
        theme_minimal() +
        theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
              axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 15),
              axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15),
              legend.background = element_rect(colour = "black", fill = "white"),
              legend.title = element_text(size = 15), legend.spacing.y = unit(.75, "cm"),
              legend.text = element_text(size = 13), legend.key.size = unit(1.5,"line"))  +
        scale_x_continuous(name = "Test Statistic") + ggtitle("Estimated Permutation Densities")
      return(list(perm_dist = data.frame("CvM" = cvm_smooth, "KS" = ks_smooth, "D1" = D1, "D2" = D2, "D3" = D3),
                  "obs_stats" = c("CvM" = cvm_obs, "KS" = ks_obs, "Diff_1" = D1_obs, "Diff_2" = D2_obs, "Diff_3" = D3_obs),
                  "Density_Plot" = gdensity,  "Predicted_Plot" = gpreds, "Permuted_Plot" = gperm))
    }
    df_preds <- data.frame("preds" = c(dfo1, dfo2), "yrs" <- c(rep("0809", 166), rep("1013", 166)))
    gpreds <- ggplot(data = df_preds, aes(x = rep(200:365,2), y = preds, color = yrs, linetype = yrs)) + geom_line(lwd = 1.15) +
      scale_color_manual(name = "Training\nYears", labels = c("2008-2009","2010-2013"), values = gg_color_hue(2)) + 
      scale_linetype_manual(name = "Training\nYears", labels = c("2008-2009","2010-2013"), values = c("solid", "dashed")) + 
      coord_cartesian(xlim = c(200, 365), ylim = c(0, max.y)) +
      scale_x_continuous(name = "Day of Year", breaks = seq(200, 365, 15)) + 
      scale_y_continuous(name = "Occurrence", breaks = seq(0, max.y, .05)) + ggtitle("Observed Prediction Curves") + theme_classic()+
      theme(plot.title = element_text(size = 20, hjust = .5), axis.title = element_text(size = 17),
            axis.text = element_text(size = 17), legend.background = element_rect(colour = "black", fill = "white"),
            legend.title = element_text(size = 17), legend.spacing.y = unit(1, "cm"),
            legend.text = element_text(size = 15),legend.key.size = unit(2,"line")) 
    
    df_perm <- data.frame("preds" = c(dfp1[index]$p0809s$y, dfp2[index]$p1013s$y), 
                          "yrs" <- c(rep("0809", 166), rep("1013", 166)))
    gperm <- ggplot(data = df_perm, aes(x = rep(200:365,2), y = preds, color = yrs, linetype = yrs)) + geom_line(lwd = 1.15) +
      scale_color_manual(name = "Training\nYears", labels = c("2008-2009","2010-2013"), values = gg_color_hue(2)) + 
      scale_linetype_manual(name = "Training\nYears", labels = c("2008-2009","2010-2013"), values = c("solid", "dashed")) + 
      coord_cartesian(xlim = c(200, 365), ylim = c(0, max.y)) +
      scale_x_continuous(name = "Day of Year", breaks = seq(200, 365, 15)) + 
      scale_y_continuous(name = "Occurrence", breaks = seq(0, max.y, .05)) + ggtitle("Permuted Prediction Curves") + theme_classic()+
      theme(plot.title = element_text(size = 20, hjust = .5), axis.title = element_text(size = 17),
            axis.text = element_text(size = 17), legend.background = element_rect(colour = "black", fill = "white"),
            legend.title = element_text(size = 17), legend.spacing.y = unit(1, "cm"),
            legend.text = element_text(size = 15),legend.key.size = unit(2,"line")) 
    return(list("Predicted_Plot" = gpreds, "Permuted_Plot" = gperm))
  }
}
library(reshape2)
load("C:/Users/drain/OneDrive/eBird project/p1013.RDA")
load("C:/Users/drain/OneDrive/eBird project/p0809.RDA")

# Testing stuff
set.seed(100)
stat_5 <- smoothed_stats(p0809, p1013, preds_0809, preds_1013, bw = 7.5, sample(1:1000,1), max.y = .35)
mean(stat_5$obs_stats["CvM"]< stat_5$perm_dist$CvM)
mean(stat_5$obs_stats["KS"] < stat_5$perm_dist$KS)
mean(stat_5$obs_stats["Diff_1"] < stat_5$perm_dist$D1)
mean(stat_5$obs_stats["Diff_2"]> stat_5$perm_dist$D2)
mean(stat_5$obs_stats["Diff_3"] > stat_5$perm_dist$D3)

stat_5$Density_Plot
stat_5$Predicted_Plot
stat_5$Permuted_Plot


## Figure 7 (a,b)
## Max Temp only
set.seed(10)
stat_unsmooth <- smoothed_stats(p0809, p1013, preds_0809, preds_1013, 7.5, sample(1:1000,1), unsmooth = TRUE, max.y = .35)
stat_unsmooth$Predicted_Plot
stat_unsmooth$Permuted_Plot
stat_unsmooth$Density_Plot
stat_unsmooth$D1to3Density
p_vals <- data.frame("KS" = mean(stat_unsmooth$obs_stats[2] < stat_unsmooth$perm_dist$KS), 
                     "CvM" = (mean(stat_unsmooth$obs_stats[1] < stat_unsmooth$perm_dist$CvM)),
                     "D1" = (mean(stat_unsmooth$obs_stats[3] < stat_unsmooth$perm_dist$D1)),
                     "D2" = (mean(stat_unsmooth$obs_stats[4] > stat_unsmooth$perm_dist$D2)),
                     "D3" = (mean(abs(stat_unsmooth$obs_stats[5]) < abs(stat_unsmooth$perm_dist$D3))))
stat_unsmooth$Predicted_Plot
stat_unsmooth$Permuted_Plot
print(p_vals)

## Figure 7 (c,d)
## DoY only
set.seed(10)
p0809s<- lapply(testpts, predict, object = observed_RF_dfs$rf.0809)
p1013s <- lapply(testpts, predict, object = observed_RF_dfs$rf.1013)
p0809_npermdfs <- unlist(lapply(p0809s, FUN = mean))
p1013_npermdfs <- unlist(lapply(p1013s, FUN = mean))

stat_unsmooth_dfs <- smoothed_stats(p0809_npermdfs, p1013_npermdfs, p0809.dfs, p1013.dfs, 
                                    7.5, sample(1:1000,1), unsmooth = TRUE)
p_vals_dfs <- data.frame("KS" = mean(stat_unsmooth_dfs$obs_stats[2] < stat_unsmooth_dfs$perm_dist$KS), 
                         "CvM" = (mean(stat_unsmooth_dfs$obs_stats[1] < stat_unsmooth_dfs$perm_dist$CvM)),
                         "D1" = (mean(stat_unsmooth_dfs$obs_stats[3] < stat_unsmooth_dfs$perm_dist$D1)),
                         "D2" = (mean(stat_unsmooth_dfs$obs_stats[4] > stat_unsmooth_dfs$perm_dist$D2)),
                         "D3" = (mean(abs(stat_unsmooth_dfs$obs_stats[5]) < abs(stat_unsmooth_dfs$perm_dist$D3))))
stat_unsmooth_dfs$Density_Plot
stat_unsmooth_dfs$Predicted_Plot
stat_unsmooth_dfs$Permuted_Plot
stat_smooth_dfs <- smoothed_stats(p0809_npermdfs, p1013_npermdfs, p0809.dfs, p1013.dfs, 
                                  7.5, sample(1:1000,1))
stat_smooth_dfs$Predicted_Plot
stat_smooth_dfs$Permuted_Plot


## FIGURE 6 (a,b)
### Full Model

p0809s.full<- lapply(testpts, predict, object = observed_RF_full$rf.0809)
p1013s.full <- lapply(testpts, predict, object = observed_RF_full$rf.1013)
p0809_nperm.full <- unlist(lapply(p0809s.full, FUN = mean))
p1013_nperm.full <- unlist(lapply(p1013s.full, FUN = mean))

set.seed(10)
stat_unsmooth_full <- smoothed_stats(p0809_nperm.full, p1013_nperm.full, p0809.full, p1013.full, 
                                     7.5, sample(1:1000,1), unsmooth = TRUE)
p_vals_full <- data.frame("KS" = mean(stat_unsmooth_full$obs_stats[2] < stat_unsmooth_full$perm_dist$KS), 
                          "CvM" = (mean(stat_unsmooth_full$obs_stats[1] < stat_unsmooth_full$perm_dist$CvM)),
                          "D1" = (mean(stat_unsmooth_full$obs_stats[3] < stat_unsmooth_full$perm_dist$D1)),
                          "D2" = (mean(stat_unsmooth_full$obs_stats[4] > stat_unsmooth_full$perm_dist$D2)),
                          "D3" = (mean(abs(stat_unsmooth_full$obs_stats[5]) < abs(stat_unsmooth_full$perm_dist$D3))))
stat_unsmooth_full$Density_Plot
stat_unsmooth_full$Predicted_Plot
stat_unsmooth_full$Permuted_Plot
stat_smooth_full <- smoothed_stats(p0809_nperm.full, p1013_nperm.full, p0809.full, p1013.full,
                                   7.5, sample(1:1000,1))
stat_smooth_full$Predicted_Plot
stat_smooth_full$Permuted_Plot

## Creating Table 2
library(xtable)
xp1 <- xtable(p_vals, digits = 3)
print(xp1, include.rownames = F)

xp2 <- xtable(p_vals_dfs, digits = 4)
print(xp2, include.rownames = F)

xp3 <- xtable(p_vals_full, digits = 4)
print(xp3, include.rownames = F)

####################################################################
## V LOCAL TESTING
####################################################################
## Part (I) - spatially stratifying points

setwd("C:/Users/drain/OneDrive/eBird project/GIS_Stuff")
load("ts.bcr30fall.RData")
library("rgdal")
library("sp")
library("rworldmap")
library("randomForest")
library("maps")
library(raster)
library(FNN)
## Spatial Clipping
FWS <- readOGR("FWSApproved.shp") #Load the WCR areas
crs1 <- proj4string(FWS) #Retrieve the Projection
llong <- cbind(ts.bcr30fall$lon, ts.bcr30fall$lat)
tspoints <- SpatialPointsDataFrame(llong, data = ts.bcr30fall, proj4string = CRS(crs1)) #Declare points as SpatialPoints class
g1 <- tspoints[FWS,]

### Adding in Average Max Temp
load("C:/Users/drain/OneDrive/eBird project/SRD_tmax_1980-2007_200-365.RData", verbose = T) ## dat
load("C:/Users/drain/OneDrive/eBird project/surfing_eBird.srd.3km.locations.RData")
loc_NE <- locations[locations$lon > -78 & locations$lat > 36.5 & locations$lat < 44,]
dMtemp <- data.frame(dat[loc_NE$LID,])
dMlist <- list(); for(i in 1:165){
  dMlist[[i]] <- data.frame(DM.day <- data.frame("mtemp" = dMtemp[,i], 
                                                 "lon" = loc_NE$lon, "lat" = loc_NE$lat))
  coordinates(dMlist[[i]])=~lon+lat
  proj4string(dMlist[[i]])=CRS(crs1) # set it to lat-long
  dMlist[[i]] = spTransform(dMlist[[i]],CRS(crs1))
  gridded(dMlist[[i]]) = TRUE
  dMlist[[i]] = raster(dMlist[[i]])
  projection(dMlist[[i]]) = CRS(crs1)
  ts.day <- ts.bcr30fall %>% mutate(DOF = round(365*dfs) - 200) %>% filter(DOF == i)
  llong <- cbind(ts.day$lon, ts.day$lat)
  extract1 <- extract(x = dMlist[[i]], y =llong)
  dMlist[[i]] <-  data.frame(ts.day, "MMtemp" = extract1)
}
ts.reco <- ldply(dMlist, data.frame)
ts.bcfa <- ts.reco[complete.cases(ts.reco),] %>% mutate(t_anomaly = temp.max - MMtemp)
save(ts.bcfa, file = "ts.bcrfall30_anomaly_raster.RDA")



# Figure 8 Code -----------------------------------------------------------

load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts.bcrfall30_anomaly_raster.RDA", verbose = T)
## Smoothing
yrly_anom <- ts.bcfa %>%  dplyr::group_by(year, DOF) %>% dplyr::summarise(tanomday = mean(t_anomaly)) %>% 
  mutate(yearIND = ifelse(year %in% 2008:2009, "08-09", "10-13")) %>% 
  group_by(yearIND, DOF) %>% dplyr::summarise(tanomday2 = mean(tanomday)) %>% 
  mutate(smoothed_anom = ksmooth(x=DOF, y = tanomday2, bandwidth = 18.5, kernel = "normal")$y)

## Smoothing overall    
ovr_anom <- ts.bcfa %>%  dplyr::group_by(DOF) %>% dplyr::summarise(tanomday = mean(t_anomaly)) %>% 
  mutate(smoothed_anom = ksmooth(x=DOF, y = tanomday, bandwidth = 18.5, kernel = "normal")$y)

## Plotting
g.anom <- ggplot(data = yrly_anom, aes( x= 200+DOF, y = smoothed_anom)) + 
  geom_line(aes(col = yearIND, linetype = yearIND), size = 1.1) +geom_hline(aes(yintercept = 0)) + theme_classic() + 
  geom_line(data = ovr_anom, aes(x = 200 + DOF, y = smoothed_anom), size = 1, col = 'gray', alpha = .75) +
  scale_x_continuous(name = "Day of Year", breaks =seq(200, 365, 15), labels = seq(200, 365, 15)) +
  scale_y_continuous(name = "Max Temp Anomaly, Degrees Celsius", breaks = seq(-8, 5, 1), labels = seq(-8, 5, 1)) + 
  scale_color_manual(name = "Year\nGrouping\n", labels = c("2008-\n2009\n", "2010-\n2013"), values = gg_color_hue(2)) + 
  scale_linetype_manual(name = "Year\nGrouping\n", labels = c("2008-\n2009\n", "2010-\n2013"), values = c("solid", "dashed")) + 
  ggtitle("Max Temp Anomaly, by year") +
  theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
        axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.title = element_text(size = 13), legend.spacing.y = unit(1, "cm"),
        legend.text = element_text(size = 12), legend.key.size = unit(2.5,"line")) 

g.anom


###################################
# Spatial Clipping
###################################

load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts.bcrfall30_anomaly_raster.RDA", verbose = T)
FWS <- readOGR("FWSApproved.shp") #Load the WCR areas
crs1 <- proj4string(FWS) #Retrieve the Projection
llong <- cbind(ts.bcfa$lon, ts.bcfa$lat)
tspoints <- SpatialPointsDataFrame(llong, data = ts.bcfa, proj4string = CRS(crs1)) #Declare points as SpatialPoints class
g1 <- tspoints[FWS,]

cll <- cbind(c(-72.597656), c(41.640078))
cp <- SpatialPoints(cll, proj4string = CRS(crs1))

## Centers of Spatial "Bubbles"
cll1 <- c(-71.213, 42.309);cll2 <- c(-72.619, 41.475);cll3 <- c(-74.113, 40.897)
cll4 <- c(-74.860, 40.044);cll5 <- c(-76.223, 39.368);cll6 <- c(-76.311, 38.479)
tLoc <- list(cll1,cll2,cll3,cll4,cll5,cll6)
tLoc <- lapply(tLoc, data.frame)
tps <- list()
t_bubbles <- list()
for(i in 1:6){tps[[i]] <- SpatialPoints(t(tLoc[[i]]), proj4string = CRS(crs1));
t_bubbles[[i]] <- gBuffer(tps[[i]], width = .3)}


ts_testing <- list()
ts_training <- g1
ts_training_loc <- list()
set.seed(36116)
for(i in 1:6){
  # Stratifying
  ts_spatial <- tspoints[t_bubbles[[i]],]
  print(dim(ts_spatial), quote = FALSE)
  ts_testing[[i]] <- ts_spatial[sample(1:dim(ts_spatial)[1], 25),]
  ts_training_loc[[i]] <- ts_spatial[-ts_testing[[i]]$FID,]
  ts_training <- ts_training[-ts_testing[[i]]$FID,]
}

ts_testinglist <- lapply(ts_testing, data.frame)
ts_training <- data.frame(ts_training)

save(ts_training, file = "ts_training_anom_ndate_raster_25pts.RDA")
save(ts_testinglist, file = "ts_testpts_anom_ndate_raster_25pts.RDA")



# Figure 9 Code -----------------------------------------------------------

library(ggmap)
library(plyr)
library(dplyr)
library(RColorBrewer)
gmap1 <- get_googlemap(center = c(lon = -73.5, lat = 40), zoom = 7)
FWS2 <- fortify(FWS)
tspts.gg <- ldply(lapply(ts_testing, data.frame), data.frame) %>% mutate(ID = as.character(sort(rep(1:6, 25))))
# tspts.gg (geom_point, mapping = aes(x = lon, y = lat, col = ID),
# alpha = .5, size = 0.5)
#writeOGR_apply <- function(poly, ...){lyr = as.character(sample(1000,1))
#writeOGR(poly, dsn = "d", layer = lyr, driver = "ESRI Shapefile", ...)}


#tbubbles.gg <- lapply(t_bubbles, writeOGR_apply)
gmap2 <- ggmap(gmap1) + 
  geom_polygon(aes(x = long, y = lat, group = group), data = FWS2, fill = 'forestgreen', alpha = .5, size = 0) + 
  geom_point(aes(x = lon, y =lat, fill = ID), data = tspts.gg, colour="black",pch=21, size=2.5, alpha = .65) +
  xlab("Longitude") + ylab("Latitude") + 
  scale_fill_discrete(name = "Testing\nLocation") + ggtitle("Locations of Testing Points")+
  theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 15),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.title = element_text(size = 13), legend.spacing.y = unit(1, "cm"),
        legend.text = element_text(size = 12)) 
gmap2




# Running the U-Statistic Tests -------------------------------------------

library(rpart)
library(MASS)
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts_testpts_anom_ndate.RDA")
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts_training_anom_ndate.RDA")

set.seed(36116)
tr_o1 <- ts_training %>% dplyr::select(c(I.stationary),time:umd.C16, occur, dfs, t_anomaly)
tr_p1 <- tr_o1 %>% mutate(t_anomaly = sample(t_anomaly, replace = F))

save(tr_o1, file = "CleanedTraining_ANOM.RDA")
save(tr_p1, file = "CleanedPTraining_ANOM.RDA")

f_clean <- function (df) df %>% dplyr::select(c(I.stationary),time:umd.C16, occur, dfs, t_anomaly)
ts_clean <- lapply(ts_testinglist, f_clean)

save(ts_clean, file = "CleanedTest_ANOM.RDA")

## Loading stuff
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CleanedTraining_ANOM.RDA")
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CleanedPTraining_ANOM.RDA")
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CleanedTest_ANOM.RDA")

set.seed(36116)
## Main_f runs Algorithm 5 of Mentch & Hooker, 2016 - calculates test-statistics
Main_f <- function(tr_o, tr_p ,test,verbose=TRUE,k=floor(sqrt(dim(tr_o1)[1])),nx1=50,nmc=150,
                   minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0) {
  # tr_o <- unpermuted dataset
  # tr_p <- permuted dataset
  # testvars -- index of columns in the data frame corresponding to variables being tested
  # verbose -- if TRUE, outputs progress
  # k -- size of each subsample
  # nx1 -- variance estimation parameter; number of estimates of conditional expectation
  # nmc -- variance estimation parameter; number of monte carlo samples used to estimate conditional expectation
  # minsplit -- the minimum number of observations needed in a leaf to be eligible for splitting
  # maxcompete -- number of "back-up" splits to keep track of
  # maxsurrogate -- number of surrogate splits to consider
  # usesurrogate -- if nonzero, this helps deal with missing data
  
  # Defining rpart Control Parameters
  control.sim <- rpart.control(minsplit=minsplit,maxcompete=maxcompete,
                               maxsurrogate=maxsurrogate,usesurrogate=usesurrogate)
  
  # Defining the size of the training set and ensemble
  n <- dim(tr_p)[1]
  m <- nx1*nmc
  test_d <- dim(test)[1]
  
  # Defining the reduced data:
  #train.red <- train[,-testvars]
  # test.red <- test[,-testvars]
  
  # Build the trees and estimate the parameters:
  pred.all <- matrix(0,nrow=1,ncol=test_d)
  diff.all <- matrix(0,nrow=1,ncol=test_d)
  cond.exp.full <- matrix(0,nrow=nx1,ncol=test_d)
  cond.exp.diff <- matrix(0,nrow=nx1,ncol=test_d)
  for (i in 1:nx1) {
    ind.x1 <- sample(1:n,size=1,replace=FALSE)
    pred.full <- matrix(0,nrow=nmc,ncol=test_d)
    pred.red <- matrix(0,nrow=nmc,ncol=test_d)
    pred.diff <- matrix(0,nrow=nmc,ncol=test_d)
    for (j in 1:nmc) {
      ind <- c(ind.x1,sample((1:dim(tr_o)[1])[-ind.x1],k-1,replace=FALSE))
      ss.full <- tr_o[ind,]	
      ss.red <- tr_p[ind,]
      tree.full <- rpart(occur~.,data=ss.full,control=control.sim)
      tree.red <- rpart(occur~.,data=ss.red,control=control.sim)
      pred.full[j,] <- predict(tree.full,test)
      pred.red[j,] <- predict(tree.red,test)
      pred.diff[j,] <- pred.full[j,] - pred.red[j,]
      #if (verbose) cat("nx1:  ",i,"          nmc:  ",j,"\n")
    }
    pred.all <- rbind(pred.all,pred.full)
    diff.all <- rbind(diff.all,pred.diff)
    cond.exp.full[i,] <- apply(pred.full,2,mean)
    cond.exp.diff[i,] <- apply(pred.diff,2,mean)
    if (verbose) cat("nx1:  ",i, "\n")
  }
  pred.all <- pred.all[-1,]
  diff.all <- diff.all[-1,]
  
  mean.full <- apply(pred.all,2,mean)
  mean.diff <- apply(diff.all,2,mean)
  
  zeta1.full <- apply(cond.exp.full,2,var)
  zeta1.diff <- cov(cond.exp.diff)
  
  zetak.full <- apply(pred.all,2,var)
  zetak.diff <- cov(diff.all)
  
  sd.full <- sqrt((m/n)*((k^2)/m)*zeta1.full + (1/m)*zetak.full)
  lbounds.full <- qnorm(0.025,mean=mean.full,sd=sd.full)
  ubounds.full <- qnorm(0.975,mean=mean.full,sd=sd.full)
  
  sd.diff <- sqrt((m/n)*((k^2)/m)*zeta1.diff + (1/m)*zetak.diff)
  lbounds.diff <- qnorm(0.025,mean=mean.diff,sd=sd.diff)
  ubounds.diff <- qnorm(0.975,mean=mean.diff,sd=sd.diff)
  
  cov.diff <- (m/n)*((k^2)/m)*zeta1.diff + (1/m)*zetak.diff
  
  tstat <- t(mean.diff) %*% ginv((m/n)*((k^2)/m)*zeta1.diff + (1/m)*zetak.diff) %*% mean.diff
  pval <- 1-pchisq(tstat,df=test_d)
  # lbounds -- lower bounds for the confidence intervals
  # ubounds -- upper bounds for the confidence intervals
  # pred -- predictions from the ensemble at the test points (center of the confidence intervals)
  # tstat -- test statistic from the test for significance
  # pval -- pvalue from the test for significance
  return(list("lbounds.diff"=lbounds.diff,"ubounds.diff"=ubounds.diff,
              "pred.diff"=mean.diff, "pred.full"=mean.full, "cov.diff" = cov.diff,
              "lbounds.full" = lbounds.full, "ubounds.full" = ubounds.full, "tstat"=tstat,"pval"=pval))
}


### Running the test

load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts_training_anom_ndate_raster_25pts.RDA", verbose = T)
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts_testpts_anom_ndate_raster_25pts.RDA", verbose = T)
library(dplyr)

set.seed(36116)
tr_o1 <- ts_training %>% dplyr::select(c(I.stationary),time:umd.C16, occur, dfs, t_anomaly)
tr_p1 <- tr_o1 %>% mutate(t_anomaly = sample(t_anomaly, replace = F))

save(tr_o1, file = "CleanedTraining_ANOM_25pts.RDA")
save(tr_p1, file = "CleanedPTraining_ANOM_25pts.RDA")

f_clean <- function (df) df %>% dplyr::select(c(I.stationary),time:umd.C16, occur, dfs, t_anomaly)
ts_clean <- lapply(ts_testinglist, f_clean)

save(ts_clean, file = "CleanedTest_ANOM_25pts.RDA")

test_results_anom_rast25pts <- list()
for(i in 1:6){
  test_results_anom_rast25pts[[i]] <- Main_f(tr_o = tr_o1, tr_p = tr_p1, test = ts_clean[[i]], 
                                             verbose=FALSE,k=160,nx1=250,nmc=5000,minsplit=10,
                                             maxcompete=0,maxsurrogate=0,usesurrogate=0)
  print(test_results_anom_rast25pts[[i]]$tstat)
  print(test_results_anom_rast25pts[[i]]$pval)
  print(i)
}

save(test_results_anom_rast25pts, file = "tresults_25pts.RDA")

## Creating Table 3

library(xtable)
save(test_results_anom_rast2, file ="eBird_anom_raster2_1.RDA")
tstats <- c();pvals <- c();for(i in 1:6){tstats[i]<- test_results_anom_rast25pts[[i]]$tstat; 
pvals[i] <- test_results_anom_rast25pts[[i]]$pval}
tstats_df<- data.frame("Testing Zone" = 1:6, tstats, pvals)
write.csv(tstats_df, "test_statistics_raster25pts.csv")
tstats_df$Testing.Zone <- as.character(1:6)
tstats_df$tstats <- as.character(round(tstats_df$tstats,2))
xstats <- xtable(tstats_df, digits = -3)
print(xstats, include.rownames = F)
