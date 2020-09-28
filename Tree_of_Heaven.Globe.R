require(ggplot2)
require(mapproj)
require(rgdal)
require(kuenm)
require(dismo)
require(rgbif)

#### Gather data ####
aa.gbif <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('Longitude', 'Latitude')
aa.bien <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
pts <- SpatialPoints(coords = unique(rbind(aa.gbif, aa.bien)))
d <- getData('worldclim', download=T, var='bio', res=5)
e <- stack(d[[1]], d[[11]], d[[12]], d[[19]])

#### Split data into testing and training sets ####
k <- 5; set.seed(1991)
folds <- kfold(pts, k)
test <- pts[which(folds!=k)]
train <- pts[which(folds==k)]

#### Model and Predict ####
m <- maxent(p=pts, x=d) #m2 <- maxent(p=pts, x=e)
p <- predict(m, x=d) #p2 <- predict(m2, x=e)

#### Map results ####
rpts <- rasterToPoints(p)
rdf <- as.data.frame(rpts)

ggc <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
  scale_fill_continuous('Suitability', type='viridis') + 
  coord_equal() + theme_void() + theme(legend.position=c(.5, .05), legend.direction = 'horizontal')

ggcp <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y), fill='grey80') + 
  geom_point(data=data.frame(coordinates(pts)), aes(x=Longitude, y=Latitude), col='red', cex=1.5, stroke=1) + 
  coord_equal() + theme_void() + theme(legend.position='none')

png('C:\\Users\\bjselige\\Desktop\\toh.global.png', width=1920, height=1080)
ggc
dev.off()

png('C:\\Users\\bjselige\\Desktop\\toh.global_points.png', width=1920, height=1080)
ggcp
dev.off()