


load("/Users/df36/projects/2017_Lucas_RF_CLT/data.surf/SRD_tmax_2014_200-365.RData")
str(dat)
load("/Users/df36/projects/2017_Lucas_RF_CLT/data.surf/surfing_eBird.srd.3km.locations.RData")
head(locations)


library(fields)
quilt.plot(x=locations$lon, y=locations$lat, z=dat[,50])
map("state", add = T)
summary(dat[,50])

   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  -30.0     4.5     7.5     7.9    11.5    27.5  687325
