require(ggplot2)
require(mapproj)
require(rgdal)
require(kuenm)
require(dismo)
require(rgbif)


# aa.gbif <- occ_data(scientificName = 'Ailanthus altissima', limit=99999)
# aa.gbif <- na.omit(data.frame(aa.gbif$data$decimalLongitude, aa.gbif$data$decimalLatitude))
# write.csv(aa.gbif, 'C:\\Users\\bjselige\\Desktop\\Ailanthus.GBIF.csv')
aa.gbif <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.GBIF.csv')
aa.gbif <- SpatialPoints(coords = aa.gbif)
aa.bien <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
aa.bien <- SpatialPoints(coords = aa.bien)

#b <- SpatialPoints(coords = na.omit(data.frame(a$data$decimalLongitude, a$data$decimalLatitude)))

b0 <- kfold(b, 5)
b1 <- b0[which(b0!=5)]
b2 <- b0[which(b0==5)]

# write.csv(b, 'C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\ku_sdm\\A_altissima\\Sp_joint.csv')
# write.csv(b1, 'C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\ku_sdm\\A_altissima\\Sp_train.csv')
# write.csv(b2, 'C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\ku_sdm\\A_altissima\\Sp_test.csv')


d <- getData('worldclim', download=T, var='bio', res=5)
e <- stack(d[[1]], d[[11]], d[[12]], d[[19]])
c <- maxent(p=b, x=d)
#c2 <- maxent(p=b, x=e)
p <- predict(c, x=d)
#p2 <- predict(c, x=e)


rpts <- rasterToPoints(p)
rdf <- as.data.frame(rpts)

ggc <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
  scale_fill_continuous('Suitability', type='viridis') + 
  coord_equal() + theme_void() + theme(legend.position=c(.5, .05), legend.direction = 'horizontal')

ggcp <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y), fill='grey80') + 
  geom_point(data=data.frame(coordinates(b)), aes(x=a.decimalLongitude, y=a.decimalLatitude), col='red', cex=1.5, stroke=1) + 
  #scale_fill_continuous(type='viridis') + 
  coord_equal() + theme_void() + theme(legend.position='none')

png('C:\\Users\\bjselige\\Desktop\\toh.global.png', width=1920, height=1080)
ggc
dev.off()


png('C:\\Users\\bjselige\\Desktop\\toh.global_points.png', width=1920, height=1080)
ggcp
dev.off()