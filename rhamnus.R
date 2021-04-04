#0.PREPARACIÓN####
packages<-c("tidyverse","biomod2","sp","rgdal",
            "raster","ggfortify","FactoMineR","ggfortify","factoextra","ecospat","hier.part","reshape")
sapply(packages,require,character.only=TRUE)

projectionUTM<-"+proj=utm +zone=28 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
projectiongeo<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

setwd("~/rhamnus")

canarias<-readOGR("islas.shp") 
canariasgeo<-spTransform(canarias, projectiongeo)

presence.absence.raster <- function (mask.raster,species.data,raster.label="") {
  require(raster)
  
  # set the background cells in the raster to 0
  mask.raster[!is.na(mask.raster)] <- 0
  
  #set the cells that contain points to 1
  speciesRaster <- rasterize(species.data,mask.raster,field=1)
  speciesRaster <- merge(speciesRaster,mask.raster)
  
  #label the raster
  names(speciesRaster) <- raster.label
  return(speciesRaster)
}

#1.DATOS CLIMÁTICOS####

setwd("~/Spalmensis")
bio01<-raster("bio01canarias.asc")
variables<-stack()
variables_list<-list.files(pattern="canarias.asc")
for(vl in variables_list){
  print(paste0("masking",vl))
  tempvariable<-raster(vl) %>%mask(canarias)
  variables<-stack(variables,tempvariable)
}

# bio02<-raster("bio02canarias.asc")%>%mask(canarias)
# bio03<-raster("bio03canarias.asc")%>%mask(canarias)
# bio04<-raster("bio04canarias.asc")%>%mask(canarias)
# #bio05<-raster("bio05canarias.asc")%>%mask(canarias)
# #bio06<-raster("bio06canarias.asc")%>%mask(canarias)
# bio07<-raster("bio07canarias.asc")%>%mask(canarias)
# bio08<-raster("bio08canarias.asc")%>%mask(canarias)
# bio09<-raster("bio09canarias.asc")%>%mask(canarias)
# #bio10<-raster("bio10canarias.asc")%>%mask(canarias)
# #bio11<-raster("bio11canarias.asc")%>%mask(canarias)
# bio12<-raster("bio12canarias.asc")%>%mask(canarias)
# #bio13<-raster("bio13canarias.asc")%>%mask(canarias)
# #bio14<-raster("bio14canarias.asc")%>%mask(canarias)
# bio15<-raster("bio15canarias.asc")%>%mask(canarias)
# #bio16<-raster("bio16canarias.asc")%>%mask(canarias)
# #bio17<-raster("bio17canarias.asc")%>%mask(canarias)
# bio18<-raster("bio18canarias.asc")%>%mask(canarias)
# bio19<-raster("bio19canarias.asc")%>%mask(canarias)
# # tpi<-raster("tpicanarias.asc")%>%mask(canarias)
# # slope<-raster("slopecanarias.asc")%>%mask(canarias)
# 
# variables<-stack(bio01,bio04,bio07,bio08,bio09,bio12,bio15,bio18,bio19)


#2.RHAMNUS####
setwd("~/rhamnus")
rhamnus<-read_delim("Rhamnus.csv",delim=";") %>% 
  #group_by(SPECIES) %>% 
  #arrange(y,by_group=TRUE) %>% 
  filter(y>1)
rhamnus<-rhamnus[c(3,6,7)]
species<-levels(factor(rhamnus$SPECIES))
coordinates(rhamnus)<-~x+y

proj4string(rhamnus)<-projectiongeo
rhamnusUTM<-spTransform(rhamnus,projectionUTM) %>% as.data.frame()

data<-raster::extract(x=variables,y=rhamnusUTM[c("x","y")],layer=1, nl=21,na.rm=TRUE)%>%data.frame()
data<-cbind(data,rhamnus$SPECIES)%>%na.omit()

colnames(data)[ncol(data)]<-"species"

pca<-PCA(data[1:ncol(data)-1],graph=F)

fviz_pca_ind(pca,
             label = "none", # hide individual labels
             habillage = as.factor(data$species), # color by groups
             palette = c("#00AFBB", "#E7B800","#33ffa8", "#FF5733"),
             addEllipses = TRUE # Concentration ellipses
)
species<-species[-c(3)]
for(s in species){
  especies<-filter(data,species==s) %>% na.omit()
  cor.rhamnus<-cor(especies[1:ncol(especies)-1],method="pearson")
  npred<-ecospat.npred(cor.rhamnus,th=0.75)
  print(paste(s,npred))
  rand.especies<-as.data.frame(sampleRandom(variables,nrow(especies)*2, ext=variables,na.rm=TRUE,xy=FALSE))
  rand.especies$pres<-rep(0,nrow(rand.especies))
  especies$pres<-rep(1,nrow(especies))
  especies["species"]<-NULL
  hier.data<-rbind(rand.especies,especies) %>% na.omit()
  random<-randomForest::randomForest(x=hier.data[1:length(hier.data)-1],y=hier.data$pres)
  seleccion<-randomForest::importance(random) %>% data.frame() %>% arrange(desc(IncNodePurity)) %>% 
    data.frame()
  seleccion<-rownames(seleccion)
  seleccion<-seleccion[1:12]
  hier.data<-hier.data %>% dplyr::select(seleccion,pres)
  hier<-hier.part(hier.data[,length(hier.data)],hier.data[1:length(hier.data)-1],barplot=TRUE)
  hier<-hier$I.perc
  hier<-arrange(hier,desc(ind.exp.var))
  hier<-rownames(hier)
  hier<-hier[1:npred]
  
  xy.especies<- subset(rhamnusUTM,SPECIES==s)
  xy.especies<-xy.especies[2:3]
  
  especie<-s
  
  especiesraster<-as.data.frame(presence.absence.raster(bio01,xy.especies[1:2],raster.label=s),xy=TRUE)
  especies.na<-especiesraster[,3]
  especies.na[especies.na==0]="NA"
  especies.na<-cbind(especiesraster[,1:2],especies.na)
  colnames(especies.na)<-c("x","y",s)
  
    # the name of studied species
  myRespName <- s
  # the presence/absences data for our species
  myResp <- as.numeric(especies.na[,myRespName])
  myResp[myResp == 2]  <- "NA"
  myResp<-as.numeric(myResp)
  # the XY coordinates of species data
  myRespXY <- especies.na[,c("x","y")]
  predictores<-stack()
  for(h in hier){
    setwd("~/Spalmensis")
    print(paste("Cropping",h))
    capatemporal<-variables[[h]]
    
    predictores<-stack(predictores,capatemporal)
  }
  setwd("~/Rhamnus")
  datos.all<-BIOMOD_FormatingData(resp.var=myResp,
                                  expl.var=predictores,
                                  resp.xy = myRespXY,
                                  resp.name = myRespName,
                                  eval.resp.var = NULL,
                                  eval.expl.var = NULL,
                                  eval.resp.xy = NULL,
                                  PA.nb.rep = 2,
                                  PA.nb.absences = nrow(xy.especies)*15,
                                  PA.strategy = 'random',
                                  na.rm = TRUE)
  
  rm(myRespXY,myResp)
  modelos.all <- BIOMOD_Modeling(
    datos.all,
    models = c('RF',"GLM","GBM"),
    #models.options = myBiomodOption,
    NbRunEval=3, 							#linea modificada para prueba. Original:5
    DataSplit=85,
    VarImport=1,
    prevalence=0.5,
    models.eval.meth = c("TSS"),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
  
  proj.modelos.all<-BIOMOD_Projection(modelos.all,
                                      new.env=predictores,
                                      proj.name=s,
                                      selected.models = 'all',
                                      binary.meth = 'TSS',
                                      filtered.meth ='TSS',
                                      compress = TRUE,
                                      build.clamping.mask = FALSE)
  
  
  modelos.ensemble.all<- BIOMOD_EnsembleModeling( modelos.all,
                                                  chosen.models = 'all',
                                                 # em.by = 'all',
                                                  eval.metric = 'all',
                                                  eval.metric.quality.threshold = 0.7,                                        ,
                                                  models.eval.meth = 'TSS',
                                                  prob.mean = FALSE,
                                                  prob.cv = FALSE,
                                                  prob.ci = FALSE,
                                                  prob.ci.alpha = 0.05,
                                                  prob.median = FALSE,
                                                  committee.averaging = FALSE,
                                                  prob.mean.weight = TRUE,
                                                  prob.mean.weight.decay = 'proportional',
                                                  VarImport = 1)
  
  
  ensemble.pres.all<-BIOMOD_EnsembleForecasting( modelos.ensemble.all,
                                                 projection.output = proj.modelos.all,
                                                 selected.models = 'all',
                                                 binary.meth = c('TSS'),
                                                 filtered.meth = c("TSS"),
                                                 total.consensus=TRUE,
                                                 compress = TRUE)
  
  
  
  writeRaster(ensemble.pres.all@proj@val[[1]],paste0(s,".ensemble.asc"),overwrite=TRUE)
  # map.all<-ensemble.pres.all@proj@val[[1]]
  # plot(map.all)
  # contour(map.all,nlevels=2,levels=750,add=TRUE)
  # evaluations.all<-melt.array(get_evaluations(modelos.all))
  # evaluations.all<-subset(evaluations.all,X1=="TSS" & X2=="Testing.data")
  # evaluations.all$data<-c(rep("All",nrow(evaluations.all)))
  # ggplot(evaluations.all,aes(x=X1,y=value,col=X3))+stat_boxplot(geom="errorbar")+geom_boxplot()+theme_bw()
  # 
  # ggplot(evaluations,aes(x=X1,y=value,col=X5))+stat_boxplot(geom="errorbar")+geom_boxplot()+theme_bw()
  # ggplot(evaluations,aes(x=X1,y=value,col=X4))+stat_boxplot(geom="errorbar")+geom_boxplot()+theme_bw()
  # evals.means<-aggregate(value~X3,data=evaluations,FUN=mean)
  # evals.sd<-aggregate(value~X3,data=evaluations,FUN=sd)
  # evals.min<-aggregate(value~X3,data=evaluations,FUN=min)
  # evaluations.fail<-subset(evaluations,value<0.8)
  # nrow(evaluations.fail)
  # table(evaluations.fail$X3)
  }



#3.DORYCNIUM####

dorycnium<-read_delim("dorycnium.csv",delim=";")
dorycnium<-dorycnium[c(3,5:6)] %>% na.omit()

dor.species<-levels(factor(dorycnium$SPECIES))
coordinates(dorycnium)<-~x+y

proj4string(dorycnium)<-projectiongeo
dorycniumUTM<-spTransform(dorycnium,projectionUTM) %>% as.data.frame()

data<-raster::extract(x=variables,y=dorycniumUTM[c("x","y")],layer=1, nl=21,na.rm=TRUE)%>%data.frame()
data<-cbind(data,dorycnium$SPECIES)%>%na.omit()

colnames(data)[ncol(data)]<-"species"

pca<-PCA(data[1:ncol(data)-1],graph=TRUE)

fviz_pca_ind(pca,
             label = "none", # hide individual labels
             habillage = as.factor(data$species), # color by groups
             palette = c("#00AFBB","#00AFBB", "#E7B800","#33ffa8"),
             addEllipses = TRUE # Concentration ellipses
)

ggplot(data=canarias,aes(long,lat,group=group),fill="grey")+
  geom_polygon()+
  geom_point(data=dorycniumUTM,aes(x,y,group=NULL,col=SPECIES))+theme_bw()+
  theme_bw()
