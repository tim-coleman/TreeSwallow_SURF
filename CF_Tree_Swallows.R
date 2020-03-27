############################################################################
# The purpose of this R file is to implement some causal forest algorithms
#  suggested by an anonymous reviewer for our paper
# "Statistical Inference on Tree Swallow Migrations Using Random Forests"
############################################################################


# Loading Data/Packages ---------------------------------------------------
library(grf)
library(dplyr)
load("C:/Users/drain/OneDrive/eBird project/ts.bcr30fall.RData")

set.seed(367116)
ts_design <- ts.bcr30fall %>% dplyr::select(-c(FID, year, y, date, occur))
ts_occur <- ts.bcr30fall %>% dplyr::select(occur) %>% unlist()

subsample_inds <- sample(nrow(ts_design), 2500)
ts_design <- ts_design[subsample_inds,]
ts_occur <- ts_occur[subsample_inds]
X.full <- ts_design %>% select(-c(temp.max))
temp.max <- ts_design %>% select(temp.max) %>% unlist()

rf_0 <- regression_forest(X = ts_design, Y = ts_occur, honesty = F)
cf_0 <- causal_forest(X = X.full, Y = unlist(ts_occur),
                             W = temp.max)

## Plotting the estimated treatment effect of temperature - causal forest
plot(cf_0$W.orig, cf_0$predictions)

delta <- .5

pred.plus = predict(rf_0, data.frame(X.full, temp.max = temp.max + delta))$predictions
pred.minus = predict(rf_0, data.frame(X.full, temp.max = temp.max - delta))$predictions
effect.hat.rf = (pred.plus - pred.minus)/(2*delta)
plot(temp.max, effect.hat.rf)





# Applying to test points from section 5 ----------------------------------
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CleanedTraining_ANOM.RDA")
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CleanedPTraining_ANOM.RDA")
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CleanedTest_ANOM_25pts.RDA", verbose = T)
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/ts_testpts_anom_ndate_raster_25pts.RDA", verbose = T)
load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CF_birds_predictions.RDA")

set.seed(367116)
causal_forest_birds <- causal_forest(X = tr_o1 %>% select(-c(t_anomaly, occur)), W = tr_o1[["t_anomaly"]], Y = tr_o1[["occur"]])

save(causal_forest_birds, file = "C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CF_birds.RDA")


# Estimating treatment effects at each point in the test sets
test_outs <- lapply(ts_testinglist, function(x){
  predict(causal_forest_birds, newdata = x[,which(names(x) %in% names(tr_o1 %>% select(-c(t_anomaly, occur))))])
})

save(test_outs, file = "C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CF_birds_predictions.RDA")


# Mapping the treatment effects

library(ggmap)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(rgdal)
library(sp)

load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CF_birds_predictions.RDA")
#load("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/CF_birds.RDA")


### GOOGLE API STUFF - REDACTED FOR CODE SUBMISSION
register_google("AIzaSyBT4_l4Lwd3wJZFsVLlZxm_vWzPA5ZcUaU")

### GETTING THE MAP
gmap1 <- get_googlemap(center = c(lon = -73.5, lat = 40), zoom = 7)


### LOADING WILDLIFE CONSERVATION AREAS

FWS <- readOGR("C:/Users/drain/OneDrive/eBird project/GIS_Stuff/FWSApproved.shp") #Load the WCR areas
FWS2 <- fortify(FWS)
tspts.gg <- ldply(lapply(ts_testinglist, data.frame), data.frame) %>%
  mutate(ID = as.character(sort(rep(1:6, 25))), AvgTempEffect = unlist(test_outs))
# tspts.gg (geom_point, mapping = aes(x = lon, y = lat, col = ID),
# alpha = .5, size = 0.5)
#writeOGR_apply <- function(poly, ...){lyr = as.character(sample(1000,1))
#writeOGR(poly, dsn = "d", layer = lyr, driver = "ESRI Shapefile", ...)}



# Mapping
gmap2 <- ggmap(gmap1) + geom_polygon(aes(x = long, y = lat, group = group), data = FWS2, fill = 'forestgreen', alpha = .35, size = 0) +
  geom_point(aes(x = lon, y =lat, fill = ID, size = abs(AvgTempEffect)), alpha = 0.6, data = tspts.gg, colour="black",pch=21) +
  xlab("Longitude") + ylab("Latitude") + 
  scale_size_continuous(name = expression(paste("|", hat(tau)(x), "|"))) + #scale_size_continuous(name = "Causal Forest\nEstimated Absolute\nTreatment Effect") +
  scale_fill_discrete(name = "Testing\nLocation") + ggtitle("Estimated Effect for Testing Points")+
  theme(plot.title = element_text(lineheight=2, face="bold", hjust = .5, size = 17),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 15),
        legend.background = element_rect(colour = "black", fill = "white"),
        legend.title = element_text(size = 13), legend.spacing.y = unit(0.25, "cm"),
        legend.text = element_text(size = 12))
gmap2

pdf("C:/Users/drain/Box Sync/Ebird Tree Swallow Project/CausalForestMap.pdf", width = 7, height = 6.5)
gmap2
dev.off()


# Plotting by time
# modifying factor name to be more informative
tspts.gg <- ldply(lapply(ts_testinglist, data.frame), data.frame) %>%
  mutate(ID = paste("Zone ", as.character(sort(rep(1:6, 25)), sep = "")), AvgTempEffect = unlist(test_outs))
CF_tseries <- ggplot(data = tspts.gg) + geom_line(aes(x = dfs*365, y = AvgTempEffect), alpha = 0.75, size = 1.25) +
  geom_point(aes(x = dfs*365, y = AvgTempEffect), size = 1.8, alpha = 0.8) + facet_wrap(ID~.) + 
  xlab("Day of Year") + ylab(expression( hat(tau(x)))) + ggtitle("Estimated Effect over Time") + 
  geom_hline(aes(yintercept = 0), col = 'blue1', lty = 'dashed', alpha = 0.5, size = 1.25) +  theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5), plot.title = element_text(lineheight =2, face = "bold", hjust = .5, size = 13))
plot(CF_tseries)

pdf("C:/Users/drain/Box Sync/Ebird Tree Swallow Project/CausalForestSeries.pdf", width = 7, height = 6.5)
CF_tseries
dev.off()


###### Making a histogram of random forest predictions for response to Reviewer 2
set.seed(367116)
reg_forest_0 <- regression_forest(X = tr_o1 %>% select(-c( occur)), Y = tr_o1[["occur"]], honesty = F)
reg_forest_outs <- lapply(ts_testinglist, function(x){
  predict(reg_forest_0, newdata = x[,which(names(x) %in% names(tr_o1 %>% select(-c(occur))))])
}) %>% unlist()

pdf("C:/Users/drain/Box Sync/Ebird Tree Swallow Project/Pred_Hist.pdf", width = 5, height = 4)
ggplot(data = data.frame(`Predicted Occurrence` = reg_forest_outs)) + 
  geom_histogram(aes(x = Predicted.Occurrence), bins = 25, col = 'blue', fill = 'tomato') + 
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + ylab("") + ggtitle("Predicted Occurrence, test points from Section 5") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
dev.off()


sum(reg_forest_outs > 0.5)
