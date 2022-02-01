require(biomod2) 
require(abind)
require(ade4)
require(caret)
require(checkmate)
require(dismo)
require(doParallel)
require(dplyr)
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

source('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\getVars.R')


#usa <- readOGR('C:\\Users\\bjselige\\Downloads\\us_lower_48_states.shp')
usa <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_states.shp')
usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
# usa.sel <- usa[usa$STATE_NAME%in%c('District of Columbia', 'Delaware', 'New Jersey', 'Maryland',
#                                    'West Virginia', 'Kentucky', 'Virginia', 'Ohio', 'Tennessee',
#                                    'Pennsylvania', 'New York', 'North Carolina'),]

co.us <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_counties.shp')
co.us <- spTransform(co.us, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
co.va <- co.us[co.us$STATE_NAME=='Virginia',]; co.sel <- co.va[which(co.va$NAME%in%c('Winchester', 'Frederick')),] 

# 1. load our species data
aa.PDA <- readOGR('C:\\Users\\bjselige\\Downloads\\2019_AllStates_Visual_SurveyDownload_11_04_2019.shp')
aa.PDA <- aa.PDA[which(aa.PDA$HostStatus=='Ailanthus altissima'), 0]
aa.PDA <- spTransform(aa.PDA, crs(usa))
aa.PDA <- data.frame(aa.PDA@coords); colnames(aa.PDA) <- c('Longitude', 'Latitude')
aa.BIEN <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.BIEN.csv')
aa.BIEN <- aa.BIEN[, c('Longitude', 'Latitude')]
aa.EDD <- read.csv('H:\\Shared drives\\APHIS  Projects\\PoPS\\Case Studies\\spotted_latternfly\\Host_mapping\\Ailanthus.eddmaps.csv')
aa.EDD <- aa.EDD[-which(is.na(aa.EDD$Longitude)),]
aa.EDD <- aa.EDD[which(aa.EDD$OccStatus=="Positive"), c('Longitude', 'Latitude')]
aa.pts <- SpatialPoints(unique(rbind(aa.PDA, aa.BIEN, aa.EDD)), CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
aa.pts <- crop(aa.pts, co.sel)

# load the environmental raster layers (could be any supported format by the raster package)
# Environmental variables extracted from Worldclim 
#myExpl <- raster::getData('worldclim', download=T, var='bio', res=10)
# biodir <- 'H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
# biovars <- stack(lapply(X=list.files(biodir), FUN=function(X){raster(paste(biodir, X, sep=''))}))
# myExpl <- crop(stack(biovars[[1]], biovars[[6]], biovars[[12]]), extent(co.va))
#myExpl <- stack(aggregate(myExpl, 2))

envi.30 <- getVars(30)
myExpl <- crop(envi.30, extent(co.sel))
#myExpl <- aggregate(myExpl, 8)
myExpl <- stack(raster::mask(myExpl, co.sel))
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
myBiomodData
plot(myBiomodData)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Computing the models
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c(#'CTA', 'SRE'
                                      'GLM',
                                      'GAM',
                                      'MARS',
                                      'FDA',
                                      'GBM',
                                      'RF',
                                      'ANN',
                                      #'MAXENT.Phillips',
                                      'MAXENT.Phillips.2'),
                                    models.options = myBiomodOption,
                                    NbRunEval = 1, #5
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
                                      eval.metric = c('TSS'),
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

myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj)

plot(myBiomodEF) # reduce layer names for plotting convegences

#######
rast <- myBiomodEF@proj@val
rpts <- rasterToPoints(rast)
rdf <- as.data.frame(rpts)

ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
  scale_fill_continuous(type='viridis') + coord_equal() +
  theme_void() + theme(legend.position='none')

png(paste('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\Figures\\cova.', 
          gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
    height=1080, width=1080); plot(ggsdm); dev.off()