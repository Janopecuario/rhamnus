#0.PREPARACIÓN####
packages<-c("tidyverse","biomod2","ENMTools","sp","rgdal")
sapply(packages,require,character.only=TRUE)

projectionUTM<-"+proj=utm +zone=28 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
projectiongeo<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
canarias<-readOGR("islas.shp") 
canarias<-spTransform(canarias, projectiongeo)

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

bio01<-crop(raster("bio01canarias.asc"),template)%>%mask(mask)
bio02<-crop(raster("bio02canarias.asc"),template)%>%mask(mask)
bio03<-crop(raster("bio03canarias.asc"),template)%>%mask(mask)
bio04<-crop(raster("bio04canarias.asc"),template)%>%mask(mask)
bio05<-crop(raster("bio05canarias.asc"),template)%>%mask(mask)
bio06<-crop(raster("bio06canarias.asc"),template)%>%mask(mask)
bio07<-crop(raster("bio07canarias.asc"),template)%>%mask(mask)
bio08<-crop(raster("bio08canarias.asc"),template)%>%mask(mask)
bio09<-crop(raster("bio09canarias.asc"),template)%>%mask(mask)
bio10<-crop(raster("bio10canarias.asc"),template)%>%mask(mask)
bio11<-crop(raster("bio11canarias.asc"),template)%>%mask(mask)
bio12<-crop(raster("bio12canarias.asc"),template)%>%mask(mask)
bio13<-crop(raster("bio13canarias.asc"),template)%>%mask(mask)
bio14<-crop(raster("bio14canarias.asc"),template)%>%mask(mask)
bio15<-crop(raster("bio15canarias.asc"),template)%>%mask(mask)
bio16<-crop(raster("bio16canarias.asc"),template)%>%mask(mask)
bio17<-crop(raster("bio17canarias.asc"),template)%>%mask(mask)
bio18<-crop(raster("bio18canarias.asc"),template)%>%mask(mask)
bio19<-crop(raster("bio19canarias.asc"),template)%>%mask(mask)

tpi<-raster("tpicanarias.asc")%>%crop(template)%>%mask(mask)
slope<-raster("slopecanarias.asc")%>%crop(template)%>%mask(mask)

variables<-stack(bio01,bio02,bio03,bio04,bio05,bio06,bio07,bio08,bio09,bio10,
                 bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19,tpi,slope)

data<-raster::extract(x=variables,y=xy.sambucus,layer=1, nl=21,na.rm=TRUE)%>%data.frame()
data<-cbind(data,origen)%>%na.omit()
colnames(data)<-c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09","bio10","bio11","bio12","bio13","bio14",
                  "bio15","bio16","bio17","bio18","bio19","tpi","slope","origen")
df<-data[1:21]

autoplot(prcomp(df),data=data,
         colour="origen",
         loadings.label = TRUE, loadings.label.size = 3)

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


#2.RHAMNUS####

rhamnus<-read_delim("Rhamnus.csv",delim=";") %>% group_by(SPECIES) %>% 
  arrange(y,by_group=TRUE) %>% 
  filter(y>1)
species<-levels(factor(rhamnus$SPECIES))
ggplot(canarias,aes(long,lat,group=group))+
  geom_polygon()+
  geom_point(data=rhamnus,aes(x,y,group=NULL))+
  theme_bw()+ theme(panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"))

for(s in species){
  especies<-filter(rhamnus,SPECIES==s)
  print(s)
  
}

#3.DORYCNIUM####

dorycnium$latitude <- dms2dec(dorycnium$y)

dorycnium<-read_delim("dorycnium.csv",delim=";")