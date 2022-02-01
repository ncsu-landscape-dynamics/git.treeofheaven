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

##### Gather and prep data #####
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

envi.30 <- getVars(res=30); envi.1k <- getVars(res=1000)

toh.PDA <- read.csv('H:\\Shared drives\\Data\\Table\\Regional\\Ailanthus.PDA.csv')[, c('Longitude', 'Latitude')]
toh.BIEN <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
toh.EDD <- read.csv('H:\\Shared drives\\Data\\Table\\USA\\Ailanthus.eddmaps.csv')
toh.EDD <- toh.EDD[-which(is.na(toh.EDD$Longitude)),]
toh.EDD <- toh.EDD[which(toh.EDD$OccStatus=="Positive"), c('Longitude', 'Latitude')]
toh.FIA <- read.csv('H:\\Shared drives\\Data\\Table\\Regional\\toh.fiaNE.csv')[,c(2,3)]; colnames(toh.FIA) <- c('Longitude', 'Latitude')
toh.xy <- SpatialPoints(unique(rbind(toh.PDA, toh.BIEN, toh.EDD, toh.FIA)), CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
toh.xy <- crop(toh.xy, borders)

toh.ras <- rasterize(x=toh.xy, y=envi.1k$CAN, fun='count', background=NA)
toh.rtp <- rasterToPoints(toh.ras, spatial=T)
toh.pts <- SpatialPointsDataFrame(coords=toh.rtp@coords,
                                  data=data.frame('X'=toh.rtp$layer, 
                                                  'MAT'=extract(x=envi.1k$MAT, y=toh.rtp),
                                                  'CAN'=extract(x=envi.1k$CAN, y=toh.rtp),
                                                  'RNR'=extract(x=envi.1k$RNR, y=toh.rtp),
                                                  'POP'=extract(x=envi.1k$POP, y=toh.rtp)))

k <- 5; set.seed(1991); folds <- kfold(toh.pts, k=k)
test <- toh.pts[which(folds==k),]; train <- toh.pts[which(folds!=k),]
backg <- spsample(borders, type='regular', n=length(test))

##### Model #####
m1.gam <- gam(formula=X~MAT+CAN+RNR+POP, data=train@data)
p1.gam <- predict(envi.1k, m1.gam, type='response')
eglm <- evaluate(p=test, a=backg, model=m1.gam, x=envi.1k)

p1.gam2 <- prediction(predictions=extract(p1.gam, toh.pts)[!is.na(extract(p1.gam, toh.pts))],
                      labels=as.numeric(toh.pts$X>0)[!is.na(extract(p1.gam, toh.pts))])
e1.gam <- performance(p1.gam, "auc")

m1.rf <- randomForest(formula=X~MAT+CAN+RNR+POP, data=train, na.action = na.omit)
p1.rf <- predict(envi.1k, m1.rf)
e1.rf <- evaluate(p=test, a=backg, model=m1.rf, x=envi.1k)

mglm <- glm(formula=X~MAT+CAN+RNR+POP, data=train@data)
pglm <- predict(envi.1k, mglm, type='response')
eglm <- evaluate(p=test, a=backg, model=mglm, x=envi.1k)

mgb <- gbm(formula=X~MAT+CAN+RNR+POP, data=train@data)
pgb <- predict.gbm(newdata=data.frame(values(envi.1k)), object=mgb, type='response', n.trees=100)
pgbr <- envi.1k$MAT; values(pgbr) <- pgb
pgbr <- mask(pgbr, borders)
evaluate(p=test, a=backg, model=mgb, x=envi.1k, n.trees=100)

# mlm <- glmer(formula=X~CAN+RNR+POP+(1|MAT), data=train@data)
# plm <- predict(envi.1k, mlm, type='response')