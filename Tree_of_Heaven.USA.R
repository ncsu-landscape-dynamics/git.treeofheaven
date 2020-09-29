### This code generates a host map prediction for Tree of Heaven across 
### the lower 48 United States for use in PoPS modeling
require(BIEN)
require(dismo)
require(sf)
require(rJava)
require(rgdal)
require(parallel)
require(ENMeval)
require(sdmvspecies)
options(java.parameters = "-Xmx8000m")

##### Gather and prep data #####
usa <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_states.shp')
usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

toh.BIEN <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
toh.EDD <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.eddmaps.csv')
toh.EDD <- toh.EDD[-which(is.na(toh.EDD$Longitude)),]
toh.EDD <- toh.EDD[which(toh.EDD$OccStatus=="Positive"), c('Longitude', 'Latitude')]
toh.xy <- SpatialPoints(unique(rbind(toh.BIEN, toh.EDD)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '))
toh.xy <- crop(toh.xy, usa)

biodir <- 'H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
biovars <- stack(lapply(X=list.files(biodir), FUN=function(X){raster(paste(biodir, X, sep=''))}))
rails.d <- raster('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\Rails_Roads\\Products_generated_from_Rails_Roads\\rails.distance.1km.tif')
rails.d <- resample(rails.d, biovars[[1]], method='bilinear'); rails.d <- rails.d + .5
roads.d <- raster('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\Rails_Roads\\Products_generated_from_Rails_Roads\\roads.distance.tif')
roads.d <- resample(roads.d, biovars[[1]], method='bilinear')

#### Split data into testing and training sets ####
k <- 5; set.seed(1991); folds <- kfold(toh.xy, k=k)
train <-  toh.xy[folds!=k,]; test <- toh.xy[folds==k,]
backg <- spsample(usa, n=length(test), type='regular')
backg.xy <- coordinates(backg); colnames(backg.xy) <- c('Longitude', 'Latitude')

#### Model and predict ####
m.all <- maxent(p=train, x=biovars)
m.top3.ro <- maxent(p=train, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d))
pm.top3.ro <- predict(m.top3.ro, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d))

#### Evaluation ####
slf <- raster('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\slf_data_redone_with_all_data_sources\\slf2019_infested.tif')
slf <- projectRaster(slf, crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '))

slf.hist <- hist(pm.top3.ro[slf$slf2019_infested>0])
per.out <- sum(pm.top3.ro[slf$slf2019_infested>0]<.2, na.rm=T)/sum(pm.top3.ro[slf$slf2019_infested>0]>0, na.rm=T)

#### Post-processing prediction ####
final.out <- pm.top3.ro
final.out[which(values(final.out)<.2)] <- 0
final.out <- rescale(final.out)
final.out <- final.out*100
# writeRaster(final.out, 'C:\\Users\\bjselige\\Desktop\\toh.USA.tif', overwrite=T)
# toh.usa <- raster('C:\\Users\\bjselige\\Desktop\\toh.USA.tif'); toh.usa <- aggregate(toh.usa, 4)

#### This block of code includes models which were tried but abandoned ####
# m.rora <- maxent(p=train, x=stack(roads.d, rails.d))
# m.all.ro <- maxent(p=train, x=stack(biovars, roads.d))
# m.all.rora <- maxent(p=train, x=stack(biovars, roads.d, rails.d))
# m.top2.ro <- maxent(p=train, x=stack(biovars[[1]], biovars[[19]], roads.d))
# m.top3.rora <- maxent(p=train, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d, rails.d))
# m2.top3.ro <- maxent(p=toh.xy, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d))
# m2.top3.rora <- maxent(p=toh.xy, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d, rails.d))
# m.top4.ro <- maxent(p=train, x=stack(biovars[[1]], biovars[[10]], biovars[[14]], biovars[[19]], roads.d))
# pm.all <- predict(m.all, x=biovars)
# pm.all.ro <- predict(m.all.ro, stack(biovars, roads.d))
# pm.all.rora <- predict(m.all.rora, stack(biovars, roads.d, rails.d))
# pm.top2.ro <- predict(m.top2.ro, x=stack(biovars[[1]], biovars[[19]], roads.d))
# pm.top3.rora <- predict(m.top3.rora, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d, rails.d))
# pm2.top3.ro <- predict(m2.top3.ro, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d))
# pm2.top3.rora <- predict(m2.top3.rora, x=stack(biovars[[1]], biovars[[14]], biovars[[19]], roads.d, rails.d))
# pm.top4.ro <- predict(m.top4.ro, x=stack(biovars[[1]], biovars[[10]], biovars[[14]], biovars[[19]], roads.d))
# hist(pm.top3.rora[slf$slf2019_infested>0])
# hist(pm2.top3.ro[slf$slf2019_infested>0])
# sum(pm.top3.rora[slf$slf2019_infested>0]<.2, na.rm=T)/sum(pm.top3.rora[slf$slf2019_infested>0]>0, na.rm=T)
# sum(pm2.top3.ro[slf$slf2019_infested>0]<.2, na.rm=T)/sum(pm.top3.ro[slf$slf2019_infested>0]>0, na.rm=T)
# sum(pm2.top3.rora[slf$slf2019_infested>0]<.2, na.rm=T)/sum(pm.top3.rora[slf$slf2019_infested>0]>0, na.rm=T)