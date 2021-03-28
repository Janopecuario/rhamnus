#0.PREPARACIÓN####
packages<-c("tidyverse","biomod2","ENMTools","sp","rgdal")
sapply(packages,require,character.only=TRUE)

projectionUTM<-"+proj=utm +zone=28 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
projectiongeo<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
canarias<-readOGR("islas.shp") 
canarias<-spTransform(canarias, projectiongeo)

#2.RHAMNUS ####

rhamnus<-read_delim("Rhamnus.csv",delim=";") %>% group_by(SPECIES) %>% 
  arrange(y,by_group=TRUE) %>% 
  filter(y>1)
species<-levels(factor(rhamnus$SPECIES))
ggplot(canarias,aes(long,lat,group=group))+
  geom_polygon()+
  geom_point(data=rhamnus,aes(x,y,group=NULL))+
  theme_bw()+ theme(panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid")

for(s in species){
  especies<-filter(rhamnus,SPECIES==s)
  print(s)
  
}

#3.DORYCNIUM ####

dorycnium$latitude <- dms2dec(dorycnium$y)

dorycnium<-read_delim("dorycnium.csv",delim=";")