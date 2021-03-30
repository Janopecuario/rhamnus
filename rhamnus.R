#0.PREPARACIÓN####
packages<-c("tidyverse","biomod2","ENMTools","sp","rgdal",
            "raster","ggfortify","FactoMineR","ggfortify")
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

ggplot(canariasgeo,aes(long,lat,group=group))+geom_polygon()+
  geom_point(data=rhamnus,aes(x,y,group=NULL,col="red"))
#1.DATOS CLIMÁTICOS####

setwd("~/Spalmensis")

bio01<-raster("bio01canarias.asc")%>%mask(canarias)
bio02<-raster("bio02canarias.asc")%>%mask(canarias)
bio03<-raster("bio03canarias.asc")%>%mask(canarias)
bio04<-raster("bio04canarias.asc")%>%mask(canarias)
#bio05<-raster("bio05canarias.asc")%>%mask(canarias)
#bio06<-raster("bio06canarias.asc")%>%mask(canarias)
bio07<-raster("bio07canarias.asc")%>%mask(canarias)
bio08<-raster("bio08canarias.asc")%>%mask(canarias)
bio09<-raster("bio09canarias.asc")%>%mask(canarias)
#bio10<-raster("bio10canarias.asc")%>%mask(canarias)
#bio11<-raster("bio11canarias.asc")%>%mask(canarias)
bio12<-raster("bio12canarias.asc")%>%mask(canarias)
#bio13<-raster("bio13canarias.asc")%>%mask(canarias)
#bio14<-raster("bio14canarias.asc")%>%mask(canarias)
bio15<-raster("bio15canarias.asc")%>%mask(canarias)
#bio16<-raster("bio16canarias.asc")%>%mask(canarias)
#bio17<-raster("bio17canarias.asc")%>%mask(canarias)
bio18<-raster("bio18canarias.asc")%>%mask(canarias)
bio19<-raster("bio19canarias.asc")%>%mask(canarias)

tpi<-raster("tpicanarias.asc")%>%mask(canarias)
slope<-raster("slopecanarias.asc")%>%mask(canarias)

variables<-stack(bio01,bio02,bio03,bio04,bio07,bio08,bio09,bio12,bio15,bio18,bio19,tpi,slope)


#2.RHAMNUS####

rhamnus<-read_delim("Rhamnus.csv",delim=";") %>% 
  #group_by(SPECIES) %>% 
  #arrange(y,by_group=TRUE) %>% 
  filter(y>1)
species<-levels(factor(rhamnus$SPECIES))
ggplot(canarias,aes(long,lat,group=group))+
  geom_polygon()+
  geom_point(data=rhamnus,aes(x,y,group=NULL))+
  theme_bw()+ theme(panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"))
coordinates(rhamnus)<-~x+y

proj4string(rhamnus)<-projectiongeo
rhamnusUTM<-spTransform(rhamnus,projectionUTM) %>% as.data.frame()

data<-raster::extract(x=variables,y=rhamnusUTM[c("x","y")],layer=1, nl=21,na.rm=TRUE)%>%data.frame()
data<-cbind(data,rhamnus$SPECIES)%>%na.omit()

colnames(data)[ncol(data)]<-"species"
df<-data[1:ncol(data)-1]

autoplot(prcomp(df),data=data,
         colour="species",
         loadings.label = TRUE, loadings.label.size = 3)


pca<-PCA(data[1:ncol(data)-1],graph=F)

fviz_pca_ind(pca,
             label = "none", # hide individual labels
             habillage = as.factor(data$species), # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE # Concentration ellipses
)
sambucus.natural<-subset(sambucus,Origen=="Natural")
head(sambucus.natural)
xy.sambucus.natural<-cbind(sambucus.natural$x,sambucus.natural$y)
data.natural<-as.data.frame(raster::extract(x=variables,y=xy.sambucus.natural,layer=1, nl=21,na.rm=TRUE))
data.natural$pres<-rep(1,nrow(data.natural))

cor.natural<-cor(data.natural,method="pearson")
ecospat.npred(cor.natural,th=0.75)
rand.natural.0<-as.data.frame(sampleRandom(variables,ncol(data.natural), ext=variables,na.rm=TRUE,xy=FALSE))
rand.natural.0$pres<-rep(0,nrow(rand.natural.0))
autoplot(prcomp(data.natural[,1:21]),data=data.natural,loadings.label=TRUE)
hier.data<-rbind(rand.natural.0,data.natural)
colnames(hier.data)<-c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09","bio10","bio11","bio12","bio13","bio14",
                       "bio15","bio16","bio17","bio18","bio19","tpi","slope","pres")
hier.data$bio06<-NULL
hier.data$bio05<-NULL
hier.data$bio13<-NULL
hier.data$bio14<-NULL
hier.data$bio17<-NULL
hier.data$bio04<-NULL
hier.data$bio16<-NULL
hier.data$bio07<-NULL
hier.data$bio19<-NULL
hier.data$tpi<-NULL
hier<-hier.part(hier.data[,length(hier.data)],hier.data[1:length(hier.data)-1],barplot=TRUE)

for(s in species){
  especies<-filter(rhamnus,SPECIES==s)
  print(s)
  
}

#3.DORYCNIUM####

dorycnium$latitude <- dms2dec(dorycnium$y)

dorycnium<-read_delim("dorycnium.csv",delim=";")