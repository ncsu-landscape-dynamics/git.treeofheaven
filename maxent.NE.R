require(BIEN)
require(dismo)
require(sf)
require(rJava)
require(rgdal)
require(parallel)
require(raster)
require(ENMeval)
require(ROCR)

source('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\getVars.R')

usa <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_states.shp')
usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
borders <- usa[usa$STATE_NAME%in%c('District of Columbia', 'Delaware', 'New Jersey', 'Maryland', 
                                    'West Virginia', 'Kentucky', 'Virginia', 'Ohio', 'Tennessee',
                                    'Pennsylvania', 'New York', 'North Carolina'),]
bbox <- as(extent(borders), 'SpatialPolygons'); crs(bbox) <- crs(usa)

toh.PDA <- readOGR('C:\\Users\\bjselige\\Downloads\\2019_AllStates_Visual_SurveyDownload_11_04_2019.shp')
toh.PDA <- toh.PDA[which(toh.PDA$HostStatus=='Ailanthus altissima'), 0]
toh.PDA <- spTransform(toh.PDA, crs(usa))
toh.PDA <- data.frame(toh.PDA@coords); colnames(toh.PDA) <- c('Longitude', 'Latitude')
toh.BIEN <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.BIEN.csv')
toh.BIEN <- toh.BIEN[, c('Longitude', 'Latitude')]
toh.EDD <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.eddmaps.csv')
toh.EDD <- toh.EDD[-which(is.na(toh.EDD$Longitude)),]
toh.EDD <- toh.EDD[which(toh.EDD$OccStatus=="Positive"), c('Longitude', 'Latitude')]
toh.all <- SpatialPoints(unique(rbind(toh.PDA, toh.BIEN, toh.EDD)), CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
toh.all <- crop(toh.all, borders)

# # toh.PDA <- readOGR('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\raw_data\\2019_AllStates_Visual_SurveyDownload_11_04_2019.shp')
# toh.PDA <- readOGR('C:\\Users\\bjselige\\Downloads\\2019_AllStates_Visual_SurveyDownload_11_04_2019.shp')
# toh.PDA <- toh.PDA[which(toh.PDA$HostStatus=='Ailanthus altissima'), 0]
# toh.PDA <- spTransform(toh.PDA, crs(usa))
# toh.BIEN <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.BIEN.csv')
# toh.BIEN <- toh.BIEN[, c('Longitude', 'Latitude')]
# toh.BIEN[, 3] <- 'BIEN'; colnames(toh.BIEN)[3] <- 'Source'
# toh.EDD <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.eddmaps.csv')
# toh.EDD <- toh.EDD[-which(is.na(toh.EDD$Longitude)),]
# toh.EDD <- toh.EDD[which(toh.EDD$OccStatus=="Positive"), c('Longitude', 'Latitude')]
# toh.EDD[, 3] <- 'EDD'; colnames(toh.EDD)[3] <- 'Source'
# #toh.BIENEDD <- unique(rbind(toh.BIEN, toh.EDD))
# toh.xy1 <- SpatialPointsDataFrame(coords=toh.BIENEDD[, c(1,2)], data=data.frame(toh.BIENEDD[, 3]), proj4string = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '))
# toh.xy1 <- crop(toh.xy1, borders)
# toh.xy2 <- crop(as(toh.PDA, 'SpatialPoints'), borders)
# toh.xy3 <- SpatialPoints(coords=unique(coordinates(toh.xy1+toh.xy2)), proj4string = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '))
# toh.PDA <- data.frame(toh.PDA@coords); colnames(toh.PDA) <- c('Longitude', 'Latitude')
# toh.PDA[, 3] <- 'PDA'; colnames(toh.PDA)[3] <- 'Source'
# toh.BIEN <- SpatialPointsDataFrame(toh.BIEN[, c(1,2)], data=data.frame(toh.BIEN[, 3]),
#                                    proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '))
# toh.EDD <- SpatialPointsDataFrame(toh.EDD[, c(1,2)], data=data.frame(toh.EDD[, 3]),
#                                    proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '))
# writeOGR(toh.PDA, 'C:\\Users\\bjselige\\Desktop\\toh.PDA.shp', driver='ESRI Shapefile', layer=1)
# writeOGR(toh.BIEN, 'C:\\Users\\bjselige\\Desktop\\toh.BIEN.shp', driver='ESRI Shapefile', layer=1)
# writeOGR(toh.EDD, 'C:\\Users\\bjselige\\Desktop\\toh.EDD.shp', driver='ESRI Shapefile', layer=1)
# writeOGR(toh.xy, 'C:\\Users\\bjselige\\Desktop\\toh.NE.shp', driver='ESRI Shapefile', layer=1)

envi.30 <- getVars(30)
envi.1k <- getVars(1000)

rnr1 <- mosaic(roads.30m, rails.30m, fun='min')

fiadf <- read.csv('H:\\Shared drives\\Data\\Table\\Regional\\toh.fiaNE.csv', stringsAsFactors = F)
fia.pts <- SpatialPoints(coords=fiadf[,c(2,3)], CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
fia.pts <- crop(fia.pts, borders)

toh.set1 <- toh.all
toh.set2 <- fia.pts
toh.set3 <- toh.set1+toh.set2

k <- 5

set.seed(1991); folds.set1 <- kfold(toh.set1, k=k)
set.seed(1991); folds.set2 <- kfold(toh.set2, k=k)
set.seed(1991); folds.set3 <- kfold(toh.set3, k=k)

test.set1 <- toh.set1[which(folds.set1==k),]
test.set2 <- toh.set2[which(folds.set2==k),]
test.set3 <- toh.set3[which(folds.set3==k),]

train.set1 <- toh.set1[which(folds.set1!=k),]
train.set2 <- toh.set2[which(folds.set2!=k),]
train.set3 <- toh.set3[which(folds.set3!=k),]

backg.set1 <- spsample(borders, type='regular', n=length(test.set1))
backg.set2 <- spsample(borders, type='regular', n=length(test.set2))
backg.set3 <- spsample(borders, type='regular', n=length(test.set3))


set.seed(1991)
folds <- kfold(toh.all, k=k)
test <- toh.all[which(folds==k),]
train <- toh.all[which(folds!=k),]
backg <- spsample(borders, type='regular', n=length(test))
####### Modeling ########
i <- 1
mods.out <- list()
eval.out <- list()
pred.out <- list()
for(i in i:k){print(i)
  test <- toh.all[which(folds==i),]
  train <- toh.all[which(folds!=i),]
  m.30m.mat.4 <- maxent(p=train, x=envi.30); print('Modeled')
  e.30m.mat.4 <- evaluate(p=test, a=backg, model=m.30m.mat.4, x=envi.30); print('Evaluated')
  p.30m.mat.4 <- predict(m.30m.mat.4, x=envi.30); print('Predicted')
  mods.out[[i]] <- m.30m.mat.4
  eval.out[[i]] <- e.30m.mat.4
  pred.out[[i]] <- p.30m.mat.4
}

mean(unlist(lapply(X=c(1:5), FUN=function(X){eval.out[[X]]@auc})), na.rm=T)

#m1 <- maxent(p=train, x=stack(bio1.30m, rnr1, pop.30m, canop)); print('Modeled')
#m2 <- maxent(p=train, x=stack(bio1.30m, rnr2, pop.30m, canop)); print('Modeled')
m1 <- maxent(p=train.set1, x=envi.1k); print('Modeled')
m2 <- maxent(p=train.set2, x=envi.1k); print('Modeled')
m3 <- maxent(p=train.set3, x=envi.1k); print('Modeled')
#m4 <- maxent(p=train.set1, x=envi.30); print('Modeled')

p1 <- predict(m1, x=envi.1k)
p2 <- predict(m2, x=envi.1k)
p3 <- predict(m3, x=envi.1k)

#e1 <- evaluate(p=test, a=backg, model=m1, x=stack(bio1.30m, rnr1, pop.30m, canop)); print('Evaluated')
#e2 <- evaluate(p=test, a=backg, model=m2, x=stack(bio1.30m, rnr2, pop.30m, canop)); print('Evaluated')
e1 <- evaluate(p=test.set1, a=backg.set1, model=m1, x=envi.1k); print('Evaluated')
e2 <- evaluate(p=test.set2, a=backg.set2, model=m2, x=envi.1k); print('Evaluated')
e3 <- evaluate(p=test.set3, a=backg.set3, model=m3, x=envi.1k); print('Evaluated')



# m.all <- maxent(p=toh.xy1, x=biovars)
# p.all <- predict(m.all, x=biovars)
# e.all <- evaluate(p=toh.xy1, a=backg1, model=m.all, x=biovars)
# 
# m2.all <- maxent(p=toh.xy2, x=biovars)
# p2.all <- predict(m2.all, x=biovars)
# e2.all <- evaluate(p=toh.xy2, a=backg2, model=m2.all, x=biovars)
# 
# m3.all <- maxent(p=toh.all, x=biovars)
# p3.all <- predict(m3.all, x=biovars)
# e3.all <- evaluate(p=toh.all, a=backg3, model=m3.all, x=biovars)
# 
# summary(lm(wc2.0_bio_30s_01~wc2.0_bio_30s_05, data=data.frame(na.omit(values(biovars)))))
# summary(lm(wc2.0_bio_30s_01~wc2.0_bio_30s_06, data=data.frame(na.omit(values(biovars)))))
#
# m.all.all <- maxent(p=toh.xy1, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# p.all.all <- predict(m.all.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# e.all.all <- evaluate(p=toh.xy1, a=backg1, model=m.all.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
#   
# m2.all.all <- maxent(p=toh.xy2, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# p2.all.all <- predict(m2.all.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# e2.all.all <- evaluate(p=toh.xy2, a=backg2, model=m2.all.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# 
# m3.all.all <- maxent(p=toh.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# p3.all.all <- predict(m3.all.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# e3.all.all <- evaluate(p=toh.all, a=backg3, model=m3.all.all, x=stack(biovars, roads.1k, edge.1k, pop.1k, forest.1k, ndvi.1km, rails.1k, canop.1k))
# 
# m.mat.4 <- maxent(p=toh.xy1, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k)) #biovars[[14]], biovars[[19]],
# p.mat.4 <- predict(m.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k)) #biovars[[14]], biovars[[19]],
# e.mat.4 <- evaluate(p=toh.xy1, a=backg1, model=m.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))
# 
# m2.mat.4 <- maxent(p=toh.xy2, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))#biovars[[14]], biovars[[19]],
# p2.mat.4 <- predict(m2.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))#biovars[[14]], biovars[[19]],
# e2.mat.4 <- evaluate(p=toh.xy2, a=backg2, model=m2.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))
# 
# m3.mat.4 <- maxent(p=toh.all, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))#biovars[[14]], biovars[[19]],
# p3.mat.4 <- predict(m3.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))#biovars[[14]], biovars[[19]],
# e3.mat.4 <- evaluate(p=toh.all, a=backg3, model=m3.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, forest.1k, canop.1k))

m.mat.4 <- maxent(p=train, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, canop.1k))
p.mat.4 <- predict(m.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, canop.1k))
e.mat.4 <- evaluate(p=test, a=backg, model=m.mat.4, x=stack(biovars[[1]], roads.1k, rails.1k, pop.1k, canop.1k))

m.bio19 <- maxent(p=train, x=biovars)
p.bio19 <- predict(m.bio19, x=biovars)
e.bio19 <- evaluate(p=test, a=backg, model=m.bio19, x=biovars)

# m.250m.mat.4 <- maxent(p=train, x=stack(bio1.250m, roads.d, rails.d, popden, forest, canop.250))
# p.250m.mat.4 <- predict(m.250m.mat.4, x=stack(bio1.250m, roads.d, rails.d, popden, forest, canop.250))
# e.250m.mat.4 <- evaluate(p=test, a=backg, model=m.250m.mat.4, x=stack(bio1.250m, roads.d, rails.d, popden, forest, canop.250))

# m.30m.mat.4 <- maxent(p=train, x=stack(bio1.30m, roads.30m, rails.30m, pop.30m, forest.30m, canop))
# p.30m.mat.4 <- predict(m.30m.mat.4, x=stack(bio1.30m, roads.30m, rails.30m, pop.30m, forest.30m, canop))
# e.30m.mat.4 <- evaluate(p=test, a=backg, model=m.30m.mat.4, x=stack(bio1.30m, roads.30m, rails.30m, pop.30m, forest.30m, canop))

m2.30m.mat.4 <- maxent(p=train, x=stack(bio1.30m, roads.30m, rails.30m, pop.30m, canop))
p2.30m.mat.4 <- predict(m2.30m.mat.4, x=stack(bio1.30m, roads.30m, rails.30m, pop.30m, canop))
e2.30m.mat.4 <- evaluate(p=toh.all, a=backg, model=m2.30m.mat.4, x=stack(bio1.30m, roads.30m, rails.30m, pop.30m, canop))

#m.250 <- maxent(p=toh.all, x=stack(roads.d, rails.d, popden, forest, canop, edge))
#p.250 <- predict(m.250, x=stack(roads.d, rails.d, popden, forest, canop, edge))

#m1 <- maxent(p=toh.xy, x=stack(ndvi, phen, rails.d))
#p1 <- predict(m1, x=stack(ndvi, phen, rails.d))
#pdf <- data.frame(rasterToPoints(p1)) 



# bmp('C:\\Users\\bjselige\\Desktop\\toh_pred1.bmp', width=1440, height=1440)
# gp1
# dev.off()
# bmp('C:\\Users\\bjselige\\Desktop\\toh_pred1_wpts.bmp', width=1440, height=1440)
# gp2
# dev.off()


m2 <- maxent(p=toh.xy, x=stack(ndvi, phen, rails.d, forest))
p2 <- predict(m2, x=stack(ndvi, phen, rails.d, forest))

toh.rast <- rasterize(toh.xy, envi[[1]]*0, fun='count', background=0)


# wc_13 <- getData('worldclim', var='bio', res=0.5, download = T, lon=c(-70), lat=35)
# wc_12 <- getData('worldclim', var='bio', res=0.5, download = T, lon=c(-115), lat=35)
# wc_11 <- getData('worldclim', var='bio', res=0.5, download = T, lon=c(-125), lat=35)
# wc_23 <- getData('worldclim', var='bio', res=0.5, download = T, lon=c(-70), lat=25)
# wc_22 <- getData('worldclim', var='bio', res=0.5, download = T, lon=c(-115), lat=25)


wcm <- merge(wc_11[[1]], wc_12[[1]], wc_13[[1]], wc_22[[1]], wc_23[[1]])
wcm <- crop(wcm, extent(c(bbox[2], bbox[4], bbox[1], bbox[3]))); wcm <- mask(wcm, borders)


m1 <- maxent(p=toh.xy, x=envi)

# co.sel <- co.va[co.va$NAME%in%c('Winchester', 'Frederick'),]
# envi.va <- mask(envi, co.sel); envi.va <- crop(envi.va, extent(co.sel))

#p1 <- predict(m1, envi.va)

# if(file.exists('C:\\Users\\bjselige\\Desktop\\railsdist.tif')==F){
#   rails <- readOGR('C:\\Users\\bjselige\\Downloads\\Rails_Roads-20200117T195917Z-001\\Rails_Roads\\tl_2015_us_rails.shp')
#   rails <- as(rails, 'SpatialLines')
#   rails <- crop(rails, bbox)
#   rails <- spTransform(rails, '+proj=lcc +lat_1=38.65 +lat_2=42.88 +lat_0=37.5 +lon_0=-77.77 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
#   rast <- mask(bio1*0, borders)
#   rast <- projectRaster(rast, crs=crs('+proj=lcc +lat_1=38.65 +lat_2=42.88 +lat_0=37.5 +lon_0=-77.77 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
#   rast.p <- SpatialPoints(rasterToPoints(rast, spatial=T), proj4string = crs(rast))
#   ncores <- detectCores() - 1
#   elves <- makeCluster(ncores); clusterExport(elves, c('rast.p', 'rails')); clusterEvalQ(elves, library(rgeos))
#   dlist <- parLapply(cl=elves, X=c(1:length(rast.p)), fun=function(X){
#     x.d <- gDistance(rast.p[X], rails)
#   }); stopCluster(elves)
#   d.spdf <- SpatialPointsDataFrame(coordinates(rast.p), data=data.frame(unlist(dlist)), proj4string = crs(rast.p))
#   rails.d <- rasterize(d.spdf, rast)
#   rails.d <- projectRaster(rails.d$unlist.dlist., crs=crs(bio1))
#   writeRaster(rails.d, 'C:\\Users\\bjselige\\Desktop\\railsdist.tif')
# }
# if(file.exists('C:\\Users\\bjselige\\Desktop\\railsdist.tif')){rails.d <- raster('C:\\Users\\bjselige\\Desktop\\railsdist.tif')}

# require(raster)
# require(rgdal)
# 
# sunil <- raster('C:\\Users\\bjselige\\Downloads\\drive-download-20200427T160821Z-001\\CONUS_Maxent_Ailanthus_avg_lcc.tif')
# sunil <- projectRaster(sunil, crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
# 
# usa <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_states.shp')
# usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
# borders <- usa[usa$STATE_NAME%in%c('Delaware', 'District of Columbia', 'Maryland', 'New York', 
#                                    'Virginia', 'West Virginia', 'New Jersey', 'Pennsylvania',
#                                    'North Carolina', 'Tennessee', 'Kentucky', 'Ohio'),]
# 
# sunil <- mask(sunil, borders)
# sunil <- crop(sunil, extent(borders))
# sunil.tr <- sunil>.515