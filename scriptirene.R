#0.PREPARACIÓN####

setwd("D:/Irene_Kosteletzkya")
packages<-c("tidyverse","biomod2","sp","rgdal",
            "raster","ggfortify","FactoMineR","ggfortify","factoextra","ecospat","hier.part")
sapply(packages,require,character.only=TRUE)

kos<-read_delim("Kos_completo.csv",delim=",")
regiones<-levels(factor(kos$Region))

#leer archivos
capas<-list.files(pattern=".tif")
listacapas<-stack()
for (c in capas){
  capatemporal<-raster(c)
  listacapas<-stack(listacapas,capatemporal)
}

plot(listacapas$wc2.1_30s_bio_1)
plantilla<-listacapas$wc2.1_30s_bio_1
xy.kos<-as.data.frame(kos[5:6])
especie<-"Kos"

#1.SELECCIÓN VARIABLES####

#PCA
data<-raster::extract(x=listacapas,y=xy.kos,layer=1, nl=32,na.rm=TRUE)%>%data.frame()
data<-cbind(data,kos$Region)%>%na.omit()

colnames(data)[ncol(data)]<-"Region"
colnames(data)
radiacion<-data[20:31]
rowSums(radiacion)
data$radiacion<-rowSums(radiacion)
data[20:31]<-NULL
data<-data %>% relocate(radiacion,before=Region)
data<-data %>% relocate(Region,after=worldclim_elevation)
pca<-PCA(data[3:ncol(data)],graph=TRUE)

fviz_pca_ind(pca,
             label = "none", # hide individual labels
             habillage = as.factor(data$Region), # color by groups
             palette = c("#00AFBB", "#E7B800","#33ffa8", "#FF5733","#EA46ED"),
             addEllipses = FALSE # Concentration ellipses
)

contribuciones<-pca$var$contrib %>% data.frame()
contribuciones<-contribuciones[,1:3]
varianza<-c(0.48,0.30,0.10)

contribuciones$cont1<-contribuciones$Dim.1*varianza[1]
contribuciones$cont2<-contribuciones$Dim.2*varianza[2]
contribuciones$cont3<-contribuciones$Dim.3*varianza[3]
contribuciones[1:3]<-NULL
contribuciones$total<-rowSums(contribuciones)
contribuciones<-arrange(contribuciones,desc(total))
predictores<-rownames(contribuciones[1:8,])
#PARTICION JERARQUICA

cor.kos<-cor(data[3:ncol(data)],method="pearson")
npred<-ecospat.npred(cor.kos,th=0.75)


rand.kos<-as.data.frame(sampleRandom(listacapas,nrow(data), ext=listacapas,na.rm=TRUE,xy=FALSE))
rand.kos$pres<-rep(0,nrow(rand.kos))
rand.kos[20:31]<-NULL
rand.kos["worldclim_elevation"]<-NULL
rand.kos.pres<-data[-c(1:3)]
rand.kos.pres$pres<-rep(1,nrow(rand.kos.pres))
hier.data<-rbind(rand.kos,rand.kos.pres)
hier.data %>% select(predictores)->hier.data
hier<-hier.part(hier.data[,length(hier.data)],hier.data[1:length(hier.data)-1],barplot=TRUE)
hier<-hier$I.perc
hier<-arrange(hier,ind.exp.var)
hier<-hier[]


#2.MODELIZACIÓN####
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
predictorestif<-stack()
for (c in predictores){
  capatemporal<-raster(paste0(c,".tif"))
  predictorestif<-stack(predictorestif,capatemporal)
}

predictorestif<-aggregate(predictorestif,fac=10)
plantilla<-predictorestif$wc2.1_30s_bio_1 %>% aggregate(fac=10)
xy.kos$pres<-rep(1,nrow(xy.kos))
kosraster<-as.data.frame(presence.absence.raster(plantilla,xy.kos,raster.label=especie),xy=TRUE)

kos.na<-kosraster[,3]
kos.na[kos.na==0]="NA"
kos.na<-cbind(kosraster[,1:2],kos.na)
colnames(kos.na)<-c("x","y","Kosteletzkya")

# the name of studied species
myRespName <- "Kosteletzkya"
# the presence/absences data for our species
myResp <- as.numeric(kos.na[,myRespName])
myResp[myResp == 2]  <- "NA"
myResp<-as.numeric(myResp)
# the XY coordinates of species data
myRespXY <- kos.na[,c("x","y")]

datos.all<-BIOMOD_FormatingData(resp.var=myResp,
                                expl.var=predictorestif,
                                resp.xy = myRespXY,
                                resp.name = myRespName,
                                eval.resp.var = NULL,
                                eval.expl.var = NULL,
                                eval.resp.xy = NULL,
                                PA.nb.rep = 1,
                                PA.nb.absences = nrow(xy.kos)*12,
                                PA.strategy = 'random',
                                na.rm = TRUE)

rm(myRespXY,myResp)
modelos.all <- BIOMOD_Modeling(
  datos.all,
  models = c('RF',"GLM","MARS","GBM","GAM"),
  #models.options = myBiomodOption,
  NbRunEval=1, 							#linea modificada para prueba. Original:5
  DataSplit=85,
  VarImport=0,
  prevalence=0.5,
  models.eval.meth = c("TSS"),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

proj.modelos.all<-BIOMOD_Projection(modelos.all,
                                    new.env=predictores,
                                    proj.name='Sambucus palmensis intro',
                                    selected.models = 'all',
                                    binary.meth = 'TSS',
                                    filtered.meth ='TSS',
                                    compress = TRUE,
                                    build.clamping.mask = FALSE)


modelos.ensemble.all<- BIOMOD_EnsembleModeling( modelos.all,
                                                chosen.models = 'all',
                                                em.by = 'all',
                                                eval.metric = 'all',
                                                eval.metric.quality.threshold = 0.8,                                        ,
                                                models.eval.meth = 'TSS',
                                                prob.mean = FALSE,
                                                prob.cv = FALSE,
                                                prob.ci = FALSE,
                                                prob.ci.alpha = 0.05,
                                                prob.median = FALSE,
                                                committee.averaging = FALSE,
                                                prob.mean.weight = TRUE,
                                                prob.mean.weight.decay = 'proportional',
                                                VarImport = 0)


ensemble.pres.all<-BIOMOD_EnsembleForecasting( modelos.ensemble.all,
                                               projection.output = proj.modelos.all,
                                               selected.models = 'all',
                                               binary.meth = c('TSS'),
                                               filtered.meth = c("TSS"),
                                               total.consensus=TRUE,
                                               compress = TRUE)



writeRaster(ensemble.pres.all@proj@val[[1]],"kosteletzkiaTSS.all.asc",overwrite=TRUE)
map.all<-ensemble.pres.all@proj@val[[1]]
plot(map.all)
contour(map.all,nlevels=2,levels=750,add=TRUE)
evaluations.all<-melt.array(get_evaluations(modelos.all))
evaluations.all<-subset(evaluations.all,X1=="TSS" & X2=="Testing.data")
evaluations.all$data<-c(rep("All",nrow(evaluations.all)))
ggplot(evaluations.all,aes(x=X1,y=value,col=X3))+stat_boxplot(geom="errorbar")+geom_boxplot()+theme_bw()
lm.eval<-lm(value~X3+X4+X5,data=evaluations)
anova(lm.eval)

ggplot(evaluations,aes(x=X1,y=value,col=X5))+stat_boxplot(geom="errorbar")+geom_boxplot()+theme_bw()
ggplot(evaluations,aes(x=X1,y=value,col=X4))+stat_boxplot(geom="errorbar")+geom_boxplot()+theme_bw()
evals.means<-aggregate(value~X3,data=evaluations,FUN=mean)
evals.sd<-aggregate(value~X3,data=evaluations,FUN=sd)
evals.min<-aggregate(value~X3,data=evaluations,FUN=min)
evaluations.fail<-subset(evaluations,value<0.8)
nrow(evaluations.fail)
table(evaluations.fail$X3)