require(BIEN)
require(dismo)
require(sf)
require(rJava)
require(rgdal)
require(parallel)
require(raster)
require(spgwr)
require(GWmodel)
require(randomForest)
require(lme4)
require(gbm)
require(mgcv)
require(ROCR)
#source('getVars.R')

usa <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_states.shp')
usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
borders <- usa[usa$STATE_NAME%in%c('District of Columbia', 'Delaware', 'New Jersey', 'Maryland', 
                                   'West Virginia', 'Kentucky', 'Virginia', 'Ohio', 'Tennessee',
                                   'Pennsylvania', 'New York', 'North Carolina'),]
bbox <- as(extent(borders), 'SpatialPolygons'); crs(bbox) <- crs(usa)

biodir <- 'H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
biovars <- stack(lapply(X=list.files(biodir), FUN=function(X){raster(paste(biodir, X, sep=''))}))
biovars <- crop(biovars, extent(borders)); biovars <- mask(biovars, borders); names(biovars)[1] <- 'MAT'

getVars <- function(res){
  if(res==30){
    bio1.30m <- raster('H:\\Shared drives\\Data\\Raster\\Regional\\bio1_NE_30m.tif')
    canop <- raster('H:\\Shared drives\\Data\\Raster\\Regional\\canop_NE_30m.tif')
    pop.30m <- raster('H:\\Shared drives\\Data\\Raster\\Regional\\popden_NE_30m.tif')
    rnr.30m <-raster('H:\\Shared drives\\Data\\Raster\\Regional\\rnr_NE_30m.tif')
    out.stack <- stack(bio1.30m, canop, pop.30m, rnr.30m)
  }
  if(res==1000){
    canop.1k <- raster('H:\\Shared drives\\Data\\Raster\\Regional\\canop_NE_1k.tif')
    pop.1k <-  raster('H:\\Shared drives\\Data\\Raster\\Regional\\popden_NE_1k.tif')
    roads.1k <- raster('H:\\Shared drives\\Data\\Raster\\Regional\\roadsD_NE_1k.tif')
    rails.1k <- raster('H:\\Shared drives\\Data\\Raster\\Regional\\railsD_NE_1k.tif')
    rnr.1k <- mosaic(roads.1k, rails.1k, fun='min')
    out.stack <- stack(biovars[[1]], canop.1k, pop.1k, rnr.1k)
  }
  names(out.stack) <- c('MAT', 'CAN', 'POP', 'RNR')
  return(out.stack)
}

envi.30 <- getVars(res=30)
envi.1k <- getVars(res=1000)

toh.PDA <- readOGR('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\raw_data\\2019_AllStates_Visual_SurveyDownload_11_04_2019.shp')
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

test.set1 <- pts.set1[which(folds.set1==k),]
test.set2 <- pts.set2[which(folds.set2==k),]
test.set3 <- pts.set3[which(folds.set3==k),]

train.set1 <- pts.set1[which(folds.set1!=k),]
train.set2 <- pts.set2[which(folds.set1!=k),]
train.set3 <- pts.set3[which(folds.set1!=k),]

backg.set1 <- spsample(borders, type='regular', n=length(test.set1))
backg.set2 <- spsample(borders, type='regular', n=length(test.set2))
backg.set3 <- spsample(borders, type='regular', n=length(test.set3))

# raxy <- rasterize(x=toh.all, y=aggregate(biovars[[1]], 5), fun='count', background=0)
# raxy <- raxy*aggregate(biovars[[1]]*0+1, 5)

ras.set1 <- rasterize(x=toh.set1, y=envi.1k$CAN, fun='count', background=NA)#; ras.set1 <- ras.set1*(envi.1k$MAT*0+1)
rtp.set1 <- rasterToPoints(ras.set1, spatial=T)
pts.set1 <- SpatialPointsDataFrame(coords=coordinates(rtp.set1),
                                   data=data.frame('X'=rtp.set1$layer, 
                                                   'MAT'=extract(x=envi.1k$MAT, y=rtp.set1),
                                                   'CAN'=extract(x=envi.1k$CAN, y=rtp.set1),
                                                   'RNR'=extract(x=envi.1k$RNR, y=rtp.set1),
                                                   'POP'=extract(x=envi.1k$POP, y=rtp.set1)))

# mgam <- gam(formula=X~MAT+CAN+RNR+POP, data=pts.all@data)
# pgam <- predict(envi.1k, mgam, type='response')
# pgam2 <- prediction(predictions=extract(pgam, all.p), labels = as.numeric(all.p$layer>0))
# egam <- performance(pgam2, "auc")

m1.gam <- gam(formula=X~MAT+CAN+RNR+POP, data=train.set1@data)
p1.gam <- predict(envi.1k, m1.gam, type='response')

p1.gam2 <- prediction(predictions=extract(p1.gam, pts.set1)[!is.na(extract(p1.gam, pts.set1))],
                      labels=as.numeric(pts.set1$X>0)[!is.na(extract(p1.gam, pts.set1))])
e1.gam <- performance(p1.gam2, "auc")

m2.gam <- gam(formula=X~MAT+CAN+RNR+POP, data=pts.set2@data)
p2.gam <- predict(envi.1k, m2.gam, type='response')

m3.gam <- gam(formula=X~MAT+CAN+RNR+POP, data=pts.set3@data)
p3.gam <- predict(envi.1k, m3.gam, type='response')

m1.rf <- randomForest(formula=X~MAT+CAN+RNR+POP, data=pts.set1, na.action = na.omit)
p1.rf <- predict(envi.1k, rf)
e1.rf <- evaluate(p=test, a=backg, model=rf, x=envi.1k)

mglm <- glm(formula=X~MAT+CAN+RNR+POP, data=all.pts@data)
pglm <- predict(envi.1k, mglm, type='response')
eglm <- evaluate(p=test, a=backg, model=mglm, x=envi.1k)

mlm <- glmer(formula=X~CAN+RNR+POP+(1|MAT), data=all.pts@data)
plm <- predict(envi.1k, mlm, type='response')

# bw1 <- gwr.sel(formula=Z~MAT, data=raxy.pts)
# m1 <- gwr(formula=Z~MAT, data=raxy.pts, bandwidth = bw1, predictions=T)
# a <- rasterize(m1$SDF, raxy)
# bw2 <- bw.gwr(Z~MAT, data=raxy.pts)
# m2 <- gwr.basic(Z~MAT, data=raxy.pts, bw=bw2)

mgb <- gbm(formula=X~MAT+CAN+RNR+POP, data=all.pts@data)
pgb <- predict.gbm(newdata=data.frame(values(envi.1k)), object=mgb, type='response', n.trees=100)
pgbr <- envi.1k$MAT
values(pgbr) <- pgb
pgbr <- mask(pgbr, borders)
evaluate(p=test, a=backg, model=mgb, x=envi.1k, n.trees=100)