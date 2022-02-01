require(ggplot2)
require(mapproj)
require(rgdal)
require(kuenm)
require(dismo)
require(rgbif)

#wrld <- readOGR('H:\\Shared drives\\Data\\Vector\\Global\\TM_world_borders.gpkg')
wrld <- readOGR('C:\\Users\\bjselige\\Documents\\TM_world_borders.gpkg')

#### Gather data ####
# aa.gbif <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('Longitude', 'Latitude')
# aa.bien <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
aa.gbif <- read.csv('C:\\Users\\bjselige\\Documents\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('Longitude', 'Latitude')
aa.bien <- read.csv('C:\\Users\\bjselige\\Documents\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]

pts <- SpatialPoints(coords = unique(rbind(aa.gbif, aa.bien)))
d <- getData('worldclim', download=T, var='bio', res=10)
e <- d; d <- stack(d[[1]], d[[11]], d[[12]], d[[19]])

#### Split data into testing and training sets ####
k <- 5; set.seed(1991)
folds <- kfold(pts, k)
test <- pts[which(folds!=k)]
train <- pts[which(folds==k)]
backg <- spsample(wrld, n=length(test), type='regular')

#### Model and Predict ####
m <- maxent(p=train, x=d); m2 <- maxent(p=pts, x=e)
p <- predict(m, x=d) #p2 <- predict(m2, x=e)
e <- evaluate(p=test, a=backg, model=m, x=d)
p.tr <- p > 0.4804556

#### Map results ####
rpts <- rasterToPoints(p)
rdf <- as.data.frame(rpts)

ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
  #geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) + 
  scale_fill_continuous(type='viridis') + 
  theme_void() + theme(legend.position='none')

ggcp <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y), fill='grey80') + 
  geom_point(data=data.frame(coordinates(pts)), aes(x=Longitude, y=Latitude), col='red', cex=1.5, stroke=1) + 
  coord_equal() + theme_void() + theme(legend.position='none')

png(paste('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\Figures\\globe.maxent.', 
          gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
    height=1080, width=2160); ggsdm; dev.off()

png(width=1920, height=1080); ggcp; dev.off()