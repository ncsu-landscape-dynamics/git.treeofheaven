require(biomod2) 
require(abind)
require(ade4)
require(caret)
require(checkmate)
require(dismo)
require(doParallel)
require(dplyr)
require(plyr)
require(earth)
require(ecospat)
require(ENMeval)
require(foreach)
require(foreign)
require(gam)
require(gbm)
require(ggplot2)
require(Hmisc)
require(lattice)
require(MASS)
require(maxnet)
require(mda)
require(mgcv)
require(methods)
require(nnet)
require(parallel)
require(PresenceAbsence)
require(pROC)
require(purrr)
require(randomForest)
require(raster)
require(rasterVis)
require(reshape)
require(rlang)
require(rpart)
require(sp)
require(stats)
require(testthat)
require(tidyr)
require(utils)
require(rgdal)
require(ggplot2)
require(mapproj)
require(rgdal)
require(kuenm)
require(dismo)
require(rgbif)
require(folderfun)
require(BIEN)
require(scales)
require(ggplot2)

path.loc <- '' #Set this to where you want any outputs to be put locally
path.web <- 'H:\\Shared drives\\Data\\'


#### Gather data ####
# aa.gbif <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('Longitude', 'Latitude')
# aa.bien <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
aa.gbif <- read.csv('C:\\Users\\bjselige\\Documents\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('longitude', 'latitude')
aa.bien <- read.csv('C:\\Users\\bjselige\\Documents\\Ailanthus.BIENv4.csv')[c(2:nrow(aa.bien)), c('longitude', 'latitude')]

aa.pts <- SpatialPoints(coords = na.omit(unique(rbind(aa.gbif, aa.bien))))

# load the environmental raster layers (could be any supported format by the raster package)
# Environmental variables extracted from Worldclim 
envi <- raster::getData(name='worldclim', download=T, var='bio', res=5)
myExpl <- stack(envi[[1]], envi[[11]], envi[[12]], envi[[19]])
#bgkg <- raster('C:\\Users\\bjselige\\Downloads\\Beck_KoppenClimate\\Beck_KG_V1_present_0p083.tif')
aa.ras <- rasterize(x=aa.pts, y=myExpl[[1]], fun='count', background=0); aa.ras <- (aa.ras*(myExpl[[1]]*0+1))>0
a2.pts <- rasterToPoints(aa.ras)
myRespName <- 'A_altissima'
myResp <- a2.pts[, 3] # the presence/absences data for our species
myResp[myResp==0] <- NA # setting 'true absences' to undefined
myRespXY <- a2.pts[, c(1,2)] # the XY coordinates of species data

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.strategy = 'random',
                                     PA.nb.absences = sum(myResp, na.rm=T))
myBiomodData; plot(myBiomodData)


# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Computing the models
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c(#'CTA', 'SRE', 'MAXENT.Phillips', 'MAXENT.Phillips.2'
                                      'GLM',
                                      'GAM',
                                      'MARS',
                                      'FDA',
                                      'GBM',
                                      'RF',
                                      'ANN'),
                                    models.options = myBiomodOption,
                                    NbRunEval = 5,
                                    DataSplit = 80,
                                    Prevalence = 0.5,
                                    VarImport = 3,
                                    models.eval.meth = 'TSS',
                                    SaveObj = TRUE,
                                    rescal.all.models = FALSE,
                                    do.full.models = FALSE,
                                    modeling.id=paste(myRespName,"FirstModeling",sep=""))

myBiomodModelEval <- get_evaluations(myBiomodModelOut) # get all models evaluation
dimnames(myBiomodModelEval) # print the dimnames of this object
myBiomodModelEval[c('TSS'),"Testing.data",,,] # print the eval scores of all selected models
vars_importance <- data.frame(get_variables_importance(myBiomodModelOut)) # print variable importances
# vars_ranked <- data.frame('SRE'=as.integer(rank(vars_importance[,1])),
#                           'GLM'=as.integer(rank(vars_importance[,2])),
#                           'GAM'=as.integer(rank(vars_importance[,3])),
#                           'MARS'=as.integer(rank(vars_importance[,4])),
#                           'FDA'=as.integer(rank(vars_importance[,5])),
#                           'CTA'=as.integer(rank(vars_importance[,6])),
#                           'GBM'=as.integer(rank(vars_importance[,7])),
#                           'RF'=as.integer(rank(vars_importance[,8])),
#                           'ANN'=as.integer(rank(vars_importance[,9])),
#                           'MAXENT.Phillips'=as.integer(rank(vars_importance[,10])))
# vars_ranked[,11] <- rowSums(vars_ranked)

# 3.2 Ensembling the models
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by='all',
                                      eval.metric = 'all', #c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = F,
                                      prob.cv = F, #don't use
                                      prob.ci = F, #prob.ci.alpha = 0.05,
                                      prob.median = F,
                                      committee.averaging = F,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional' )
get_evaluations(myBiomodEM) # get evaluation scores

### 4. projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = myExpl,
                                  proj.name = 'current',
                                  selected.models = 'all',
                                  binary.meth = 'TSS',
                                  compress = 'xz',
                                  clamping.mask = F,
                                  output.format = '.RData')
plot(myBiomodProj)

myCurrentProj <- get_predictions(myBiomodProj) # if you want to make custom plots, you can also get the projected map
myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM, projection.output = myBiomodProj)
plot(myBiomodEF) # reduce layer names for plotting convegences

##### 5. Thresholding
p.out <- myBiomodEF@proj@val[[1]]
values(p.out) <- rescale(values(p.out), from=c(min(values(p.out), na.rm=T), max(values(p.out), na.rm=T)), to=c(0,1))
p2 <- p.out

trs <- seq(0.01, .99, by=.01)
trs.metrics <- ldply(.data=c(1:length(trs)), .fun=function(X){
  p2.tr <- p2>trs[X]
  x.zero <- sum(raster::extract(x=p2.tr, y=aa.pts)==0, na.rm=T)/length(aa.pts)
  x.area <- sum(values(p2.tr), na.rm=T)/length(na.omit(p2[]))
  x.score <- 1 - x.zero - x.area
  df.out <- data.frame('trs'= trs[X], 'zero'=x.zero, 'area'=x.area, 'score'=x.score, row.names = trs[X])
  return(df.out)}, .progress = 'text')

ts.gg <- ggplot(data=trs.metrics) + 
  geom_line(aes(x=trs, y=zero), col='blue') + 
  geom_line(aes(x=trs, y=area), col='red') +
  geom_line(aes(x=trs, y=score), col='orange') +
  scale_x_continuous(name='Threshold', breaks=seq(min(trs.metrics$trs), max(trs.metrics$trs), by=.02)) + 
  scale_y_continuous(name='Percentage', limits = c(0, 1)) + 
  annotate(geom='text', x=min(trs.metrics$trs)+.05, y=.02, col='blue', label='Points Omitted') +
  annotate(geom='text', x=min(trs.metrics$trs)+.05, y=.3, col='red', label='Area') +
  theme_minimal()

tr.best <- trs.metrics$trs[which.max(trs.metrics$score)]
p.tr <- p.out>tr.best

##### Save
writeRaster(p.out, filename = paste(path.loc, 'toh.global_ensemble.5.tif', sep=''), format="GTiff", overwrite=T)
writeRaster(p.tr, filename = paste(path.loc, 'toh.global_ensemble.5tr.tif', sep=''), format="GTiff", overwrite=T)

###### Plotting
# #borders <- 
# rast <- myBiomodEF@proj@val
# rpts <- rasterToPoints(rast)
# rdf <- as.data.frame(rpts)
# ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
#   #geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) + 
#   scale_fill_continuous(type='viridis') + 
#   theme_void() + theme(legend.position='none')
# png(paste('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\Figures\\globe.', 
#           gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
#     height=1080, width=2160); plot(ggsdm); dev.off()