# RASTER CODE
install.packages("scico")
library("scico")
install.packages("ggsn")
library('ggsn')

######
# Add country outline: 
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Swiz_Shapefile")
switz = readOGR("G1L12.shp", layer = "G1L12")
# reproject spatial outline to match raster data
switz_outline_RP = spTransform(switz, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

####################################################################### Read in tif, filter, extend, stack
# read in masking object
forestmask = readOGR("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Swiz_Forest_Cover_Shapefile/Vector_Landuse_CH/VEC200_LandCover.shp")
forestmask_RP = spTransform(forestmask,crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


# JULY 2018 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2018_TIF")
Swiz_rastlist_JULY_18 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2018_TIF"), 
                                    pattern='.tif$', all.files=TRUE, full.names=FALSE)

JULY_18_List = list()
for (i in Swiz_rastlist_JULY_18[1:4]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  JULY_18_List[[i]] = raster.file
}

# set model for all extend
large_dims = extent(forestmask_RP)
extend_whole = extend(JULY_18_List[[3]], large_dims)

# extend
extended_allrasters_Swiz_JULY_18 = lapply(JULY_18_List, resample, extend_whole, method = "bilinear")
stack_Swiz_JULY_18 = stack(extended_allrasters_Swiz_JULY_18)
mean_Swiz_July_18 = mean(stack_Swiz_JULY_18, na.rm = TRUE)
#plot(mask(mean_Swiz_July_18, forestmask_RP, inverse = TRUE))

# save
writeRaster(stack_Swiz_JULY_18, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/July_18_Stack",
            bylayer=TRUE,format="GTiff")

# AUG 2018 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2018_TIF")
Swiz_rastlist_AUG_18 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2018_TIF"), 
                                   pattern='.tif$', all.files=TRUE, full.names=FALSE)

AUG_18_List = list()
for (i in Swiz_rastlist_AUG_18[1:50]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  AUG_18_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_AUG_18 = lapply(AUG_18_List, resample, extend_whole, method = "bilinear")
stack_Swiz_AUG_18 = stack(extended_allrasters_Swiz_AUG_18)
mean_Swiz_AUG_18 = mean(stack_Swiz_AUG_18, na.rm = TRUE)
#plot(mask(mean_Swiz_AUG_18, forestmask_RP, inverse = TRUE))

# save
writeRaster(stack_Swiz_AUG_18, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Aug_18_Stack",
            bylayer=TRUE,format="GTiff")

# JUNE 2019 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_June_2019_TIF")
Swiz_rastlist_JUNE_19 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_June_2019_TIF"), 
                                    pattern='.tif$', all.files=TRUE, full.names=FALSE)
JUNE_19_List = list()
for (i in Swiz_rastlist_JUNE_19[1:57]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  JUNE_19_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_JUNE_19 = lapply(JUNE_19_List, resample, extend_whole, method = "bilinear")
stack_Swiz_JUNE_19 = stack(extended_allrasters_Swiz_JUNE_19)

# save
writeRaster(stack_Swiz_JUNE_19, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/June_19_Stack",
            bylayer=TRUE,format="GTiff")

# JULY 2019 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2019_TIF")
Swiz_rastlist_JULY_19 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2019_TIF"), 
                                    pattern='.tif$', all.files=TRUE, full.names=FALSE)
JULY_19_List = list()
for (i in Swiz_rastlist_JULY_19[1:18]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  JULY_19_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_JULY_19 = lapply(JULY_19_List, resample, extend_whole, method = "bilinear")
stack_Swiz_JULY_19 = stack(extended_allrasters_Swiz_JULY_19)

# save
writeRaster(stack_Swiz_JULY_19, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/July_19_Stack",
            bylayer=TRUE,format="GTiff")

# AUG 2019 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2019_TIF")
Swiz_rastlist_AUG_19 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2019_TIF"), 
                                   pattern='.tif$', all.files=TRUE, full.names=FALSE)
AUG_19_List = list()
for (i in Swiz_rastlist_AUG_19[1:52]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  AUG_19_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_AUG_19 = lapply(AUG_19_List, resample, extend_whole, method = "bilinear")
stack_Swiz_AUG_19 = stack(extended_allrasters_Swiz_AUG_19)
# what is the scale of resultion? 

# save
writeRaster(stack_Swiz_AUG_19, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Aug_19_Stack",
            bylayer=TRUE,format="GTiff")

# JUNE 2020 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_June_2020_TIF")
Swiz_rastlist_JUNE_20 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_June_2020_TIF"), 
                                    pattern='.tif$', all.files=TRUE, full.names=FALSE)

JUNE_20_List = list()
for (i in Swiz_rastlist_JUNE_20[1:48]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  JUNE_20_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_JUNE_20 = lapply(JUNE_20_List, resample, extend_whole, method = "bilinear")
stack_Swiz_JUNE_20 = stack(extended_allrasters_Swiz_JUNE_20)

# save
writeRaster(stack_Swiz_JUNE_20, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/June_20_Stack",
            bylayer=TRUE,format="GTiff")

# JULY 2020 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2020_TIF")
Swiz_rastlist_JULY_20 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2020_TIF"), 
                                    pattern='.tif$', all.files=TRUE, full.names=FALSE)
JULY_20_List = list()
for (i in Swiz_rastlist_JULY_20[1:21]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  JULY_20_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_JULY_20 = lapply(JULY_20_List, resample, extend_whole, method = "bilinear")
stack_Swiz_JULY_20 = stack(extended_allrasters_Swiz_JULY_20)

# save
writeRaster(stack_Swiz_JULY_20, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/July_20_Stack",
            bylayer=TRUE,format="GTiff")

# AUG 2020 TIF
setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2020_TIF")
Swiz_rastlist_AUG_20 <- list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2020_TIF"), 
                                   pattern='.tif$', all.files=TRUE, full.names=TRUE)
AUG_20_List = list()
for (i in Swiz_rastlist_AUG_20[1:64]) {
  raster.file = raster(i)
  values(raster.file)[values(raster.file)>4] = NA
  raster.file = aggregate(raster.file, fact = 10)
  AUG_20_List[[i]] = raster.file
}

# extend
extended_allrasters_Swiz_AUG_20 = lapply(AUG_20_List, resample, extend_whole, method = "bilinear")
stack_Swiz_AUG_20 = stack(extended_allrasters_Swiz_AUG_20)

# save
writeRaster(stack_Swiz_AUG_20, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Aug_20_Stack",
            bylayer=TRUE,format="GTiff")

# combine summers
# All Summer 18 Stack + sd + mean
stack_Swiz_Summer_18 = stack(stack_Swiz_JULY_18,stack_Swiz_AUG_18)
stack_Swiz_Summer_18_sd = calc(stack_Swiz_Summer_18, fun = sd, na.rm = TRUE)
Swiz_Summer_18_mean = mean(stack_Swiz_Summer_18, na.rm = TRUE)

# save
writeRaster(stack_Swiz_Summer_18, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_18_Stack",
            bylayer=TRUE,format="GTiff")
writeRaster(stack_Swiz_Summer_18_sd, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_18_SD",
            bylayer=TRUE,format="GTiff")
writeRaster(Swiz_Summer_18_mean, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_18_Mean",
            bylayer=TRUE,format="GTiff")

# All Summer 19 Stack + sd + mean
stack_Swiz_Summer_19 = stack(stack_Swiz_JUNE_19, stack_Swiz_JULY_19, stack_Swiz_AUG_19)
stack_Swiz_Summer_19_sd = calc(stack_Swiz_Summer_19, fun = sd, na.rm = TRUE)
Swiz_Summer_19_mean = mean(stack_Swiz_Summer_19, na.rm = TRUE)

# save
writeRaster(stack_Swiz_Summer_19, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_19_Stack",
            bylayer=TRUE,format="GTiff")
writeRaster(stack_Swiz_Summer_19_sd, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_19_SD",
            bylayer=TRUE,format="GTiff")
writeRaster(Swiz_Summer_19_mean, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_19_Mean",
            bylayer=TRUE,format="GTiff")

# All Summer 20 Stack + sd + mean
stack_Swiz_Summer_20 = stack(stack_Swiz_JUNE_20, stack_Swiz_JULY_20, stack_Swiz_AUG_20)
stack_Swiz_Summer_20_sd = calc(stack_Swiz_Summer_20, fun = sd, na.rm = TRUE)
Swiz_Summer_20_mean = mean(stack_Swiz_Summer_20, na.rm = TRUE)

#save
writeRaster(stack_Swiz_Summer_20, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_20_Stack",
            bylayer=TRUE,format="GTiff")
writeRaster(stack_Swiz_Summer_20_sd, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_20_SD",
            bylayer=TRUE,format="GTiff")
writeRaster(Swiz_Summer_20_mean, "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Summer_20_Mean",
            bylayer=TRUE,format="GTiff")

# stack all rasters and by month/year 
# only stacks work
All_Stacks = stack(stack_Swiz_JULY_18, stack_Swiz_AUG_18,
                   stack_Swiz_JUNE_19, stack_Swiz_JULY_19,stack_Swiz_AUG_19,
                   stack_Swiz_JUNE_20, stack_Swiz_JULY_20, stack_Swiz_AUG_20)

All_Stacks_18 = stack(stack_Swiz_JULY_18, stack_Swiz_AUG_18)
Mean_All_Stacks_18 = mean(All_Stacks_18, na.rm = TRUE)
All_Stacks_19 = stack(stack_Swiz_JUNE_19, stack_Swiz_JULY_19, stack_Swiz_AUG_19)
Mean_All_Stacks_19 = mean(All_Stacks_19, na.rm = TRUE)
All_Stacks_20 = stack(stack_Swiz_JUNE_20, stack_Swiz_JULY_20, stack_Swiz_AUG_20)
Mean_All_Stacks_20 = mean(All_Stacks_20, na.rm = TRUE)

All_Stacks_June = stack(stack_Swiz_JUNE_19, stack_Swiz_JUNE_20)
Mean_All_Stacks_June = mean(All_Stacks_June, na.rm = TRUE)
All_Stacks_July = stack(stack_Swiz_JULY_18,stack_Swiz_JULY_19, stack_Swiz_JULY_20)
Mean_All_Stacks_July = mean(All_Stacks_July, na.rm = TRUE)
All_Stacks_Aug = stack(stack_Swiz_AUG_18, stack_Swiz_AUG_19, stack_Swiz_AUG_20)
Mean_All_Stacks_Aug = mean(All_Stacks_Aug, na.rm = TRUE)


# find mean/sd of all rasters
SD_All_Stacks = calc(All_Stacks, fun = sd, na.rm = TRUE)
Mean_All_Stacks = mean(All_Stacks, na.rm = TRUE)

# mask then take mean/sd so its forest only 
All_Stacks_Masked = mask(All_Stacks, forestmask_RP, inverse = TRUE)
mean(values(All_Stacks_Masked), na.rm = TRUE) #1.266063
sd(values(All_Stacks_Masked), na.rm = TRUE) #0.5813665

# 2018
stack_Swiz_Summer_18_Masked = mask(stack_Swiz_Summer_18, forestmask_RP, inverse = TRUE)
mean(values(stack_Swiz_Summer_18_Masked), na.rm = TRUE) #1.246557
sd(values(stack_Swiz_Summer_18_Masked), na.rm = TRUE) #0.4965932

# 2019
stack_Swiz_Summer_19_Masked = mask(stack_Swiz_Summer_19, forestmask_RP, inverse = TRUE)
mean(values(stack_Swiz_Summer_19_Masked), na.rm = TRUE) #1.283212
sd(values(stack_Swiz_Summer_19_Masked), na.rm = TRUE) #0.6054168

# 2020
stack_Swiz_Summer_20_Masked = mask(stack_Swiz_Summer_20, forestmask_RP, inverse = TRUE)
mean(values(stack_Swiz_Summer_20_Masked), na.rm = TRUE) #1.258108
sd(values(stack_Swiz_Summer_20_Masked), na.rm = TRUE) #0.5994255

#June
All_Stacks_June_Masked = mask(All_Stacks_June, forestmask_RP, inverse = TRUE)
mean(values(All_Stacks_June_Masked), na.rm = TRUE) #1.216963
sd(values(All_Stacks_June_Masked), na.rm = TRUE) #0.5591247

#July
All_Stacks_July_Masked = mask(All_Stacks_July, forestmask_RP, inverse = TRUE)
mean(values(All_Stacks_July_Masked), na.rm = TRUE) #1.125943
sd(values(All_Stacks_July_Masked), na.rm = TRUE) #0.537743

#AUG
All_Stacks_Aug_Masked = mask(All_Stacks_Aug, forestmask_RP, inverse = TRUE)
mean(values(All_Stacks_Aug_Masked), na.rm = TRUE) #1.343962
sd(values(All_Stacks_Aug_Masked), na.rm = TRUE) #0.5977083

# mask before plot
# Mean All Summers
Mean_All_Stacks_Masked = mask(Mean_All_Stacks, forestmask_RP, inverse = TRUE)
stack_Swiz_mean_pts = rasterToPoints(Mean_All_Stacks_Masked,spatial=TRUE)
stack_Swiz_df = data.frame(stack_Swiz_mean_pts)

install.packages('ggspatial')
library('ggspatial')

All_Summers_Mean_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+ 
  labs(fill="WUE")+
  ggtitle(expression(paste("Mean Summer WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(), legend.position="bottom") +
  geom_point(data = data.corrected, aes(x = long, y = lat), shape = 24, colour = "blue", 
             size = 1, alpha = 1,) + 
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
            dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84", 
           st.dist = 0.05, st.size = 3) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5) 

median(stack_Swiz_df$layer)

# SD all Summers
SD_All_Stacks_Masked = mask(SD_All_Stacks, forestmask_RP, inverse = TRUE)
stack_Swiz_SD_pts = rasterToPoints(SD_All_Stacks_Masked,spatial=TRUE)
stack_Swiz_df_SD = data.frame(stack_Swiz_SD_pts)

All_Summers_SD_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_SD,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1)+
  labs(fill="WUE")+
  ggtitle(expression(paste("Standard Deviation of Summer WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(), legend.position="bottom") +
  geom_point(data = data.corrected, aes(x = long, y = lat), shape = 21, colour = "black",
             size = 1, alpha = 1) + 
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84", 
           st.dist = 0.05, st.size = 3) + 
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5) 

# 2 panel original image remake
WUE_SD_All_Summers_2Panel = ggarrange(All_Summers_Mean_WUE_Plot, All_Summers_SD_WUE_Plot, 
                                      labels = c("A)", "B)"), nrow = 1, ncol = 2)

# save image
ggsave(WUE_SD_All_Summers_2Panel, plot = WUE_SD_All_Summers_2Panel, device = "pdf",
       path = "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Final_Plots",
       dpi = 300, scale = 2)

# All 18 WUE
Mean_All_Stacks_18_Masked = mask(Mean_All_Stacks_18, forestmask_RP, inverse = TRUE)
stack_Swiz_18_pts = rasterToPoints(Mean_All_Stacks_18_Masked,spatial=TRUE)
stack_Swiz_df_18 = data.frame(stack_Swiz_18_pts)

All_18_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_18,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+
  ggtitle(expression(paste("Mean 2018 WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(),
        legend.position = "none") +
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84",
           st.size = 2, border.size = 0.1) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5) 

# All 19 WUE
Mean_All_Stacks_19_Masked = mask(Mean_All_Stacks_19, forestmask_RP, inverse = TRUE)
stack_Swiz_19_pts = rasterToPoints(Mean_All_Stacks_19_Masked,spatial=TRUE)
stack_Swiz_df_19 = data.frame(stack_Swiz_19_pts)

All_19_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_19,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+
  ggtitle(expression(paste("Mean 2019 WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(),
        legend.position = "none") +
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84",
           st.size = 2, border.size = 0.1) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5) 

# All 20 WUE
Mean_All_Stacks_20_Masked = mask(Mean_All_Stacks_20, forestmask_RP, inverse = TRUE)
stack_Swiz_20_pts = rasterToPoints(Mean_All_Stacks_20_Masked,spatial=TRUE)
stack_Swiz_df_20 = data.frame(stack_Swiz_20_pts)

All_20_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_20,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+
  ggtitle(expression(paste("Mean 2020 WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(),
        legend.position = "none") +
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84",
           st.size = 2, border.size = 0.1) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5) 

#f = north2(All_20_WUE_Plot, x=.9, y=.08,  symbol = 10)

# 3 panel 18, 19, 20
WUE_18_19_20_3Panel = ggarrange(All_18_WUE_Plot, All_19_WUE_Plot, All_20_WUE_Plot,
                                nrow = 1, ncol = 3)

# All June WUE
Mean_All_Stacks_June_Masked = mask(Mean_All_Stacks_June, forestmask_RP, inverse = TRUE)
stack_Swiz_June_pts = rasterToPoints(Mean_All_Stacks_June_Masked,spatial=TRUE)
stack_Swiz_df_June = data.frame(stack_Swiz_June_pts) 

All_June_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_June,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+
  ggtitle(expression(paste("Mean June WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(),
        legend.position = "none") +
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84",
           st.size = 2, border.size = 0.1) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5)

# All July WUE
Mean_All_Stacks_July_Masked = mask(Mean_All_Stacks_July, forestmask_RP, inverse = TRUE)
stack_Swiz_July_pts = rasterToPoints(Mean_All_Stacks_July_Masked,spatial=TRUE)
stack_Swiz_df_July = data.frame(stack_Swiz_July_pts)

All_July_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_July,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+
  ggtitle(expression(paste("Mean July WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(),
        legend.position = "none") +
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84",
           st.size = 2, border.size = 0.1) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5)

# All Aug WUE
Mean_All_Stacks_Aug_Masked = mask(Mean_All_Stacks_Aug, forestmask_RP, inverse = TRUE)
stack_Swiz_Aug_pts = rasterToPoints(Mean_All_Stacks_Aug_Masked,spatial=TRUE)
stack_Swiz_df_Aug = data.frame(stack_Swiz_Aug_pts)

All_Aug_WUE_Plot = ggplot()+
  geom_raster(data = stack_Swiz_df_Aug,aes(x = x, y = y, fill = layer))+
  scale_fill_scico(palette = 'batlow',direction=-1, limits=c(0,4))+
  ggtitle(expression(paste("Mean August WUE (g C ",kg^-1," ",H[2],O,")")))+
  geom_polygon(data=switz_outline_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank(),
        legend.position = "none") +
  scalebar(location = "bottomright", 
           x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
           dist = 50, dist_unit = "km", transform = TRUE,  model = "WGS84",
           st.size = 2, border.size = 0.1) +
  north(x.min = 5.958671, x.max = 10.49951, y.min = 45.81269, y.max = 47.80561,
        symbol = 10, location = "topright") + coord_fixed(ratio = 1.5)

# 3 panel June, July, Aug
WUE_June_July_Aug_3Panel = ggarrange(All_June_WUE_Plot, All_July_WUE_Plot, All_Aug_WUE_Plot,
                                nrow = 1, ncol = 3)


### end of making raster stacks

# original figure
# 18 19 20 plots - sup - masked 
# june july aug plots - sup - masked




### read in all stacks in case of failure
# AUG 18 
AUG_18_list_factored =  list.files(path = setwd("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Raster_Stacks/Aug_18_Stack"), 
                             pattern='.tif$', all.files=TRUE, full.names=FALSE)
AUG_18_List_test = list()
for (i in AUG_18_list_factored) {
  raster.file = raster(i)
  AUG_18_List_test[[i]] = raster.file
}

# add rest of files if needed







######################################################################################################

# SITE DATA CODE

# WUE ####################

# read in csv files
WUE_June_19 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_June_2019/June-2019-ECO4WUE-001-results.csv")
WUE_June_20 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_June_2020/June-2020-ECO4WUE-001-results.csv")
WUE_July_18 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_July_2018/July-2018-ECO4WUE-001-results.csv")
WUE_July_19 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_July_2019/July-2019-ECO4WUE-001-results.csv")
WUE_July_20 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_July_2020/July-2020-ECO4WUE-001-results.csv")
WUE_Aug_18 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_Aug_2018/August-2018-ECO4WUE-001-results.csv")
WUE_Aug_19 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_Aug_2019/August-2019-ECO4WUE-001-results.csv")
WUE_Aug_20 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_Aug_2020/August-2020-ECO4WUE-001-results.csv")


# rbind
all_Swiz_WUE_data = rbind(WUE_June_19, WUE_June_20, WUE_July_18, WUE_July_19, WUE_July_20,
                          WUE_Aug_18, WUE_Aug_19,WUE_Aug_20)

colnames(all_Swiz_WUE_data)[5] = "Date_total"
colnames(all_Swiz_WUE_data)[9] = "WUEavg"

# Adjust time zone + separate date and time
all_Swiz_WUE_data$date = as_datetime(all_Swiz_WUE_data$Date_total)
all_Swiz_WUE_data$date.adj = all_Swiz_WUE_data$date + hours(1)
all_Swiz_WUE_data = separate(data=all_Swiz_WUE_data, col=date.adj, into=c("date","time"), sep =" ")

# extract "m-d" from "Y-m-d"
n_last = 5
substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))
all_Swiz_WUE_data$m_d = substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))

# Make $ID character
all_Swiz_WUE_data$ID = as.character(all_Swiz_WUE_data$ID)

# seperate into seasonal
Swiz_WUE_Data_Summer_18 = all_Swiz_WUE_data[all_Swiz_WUE_data$date >= "2018-01-01" & all_Swiz_WUE_data$date <= "2018-12-31",]
Swiz_WUE_Data_Summer_19 = all_Swiz_WUE_data[all_Swiz_WUE_data$date >= "2019-01-01" & all_Swiz_WUE_data$date <= "2019-12-31",]
Swiz_WUE_Data_Summer_20 = all_Swiz_WUE_data[all_Swiz_WUE_data$date >= "2020-01-01" & all_Swiz_WUE_data$date <= "2020-12-31",]

# Read in RF_Temp_RH csv
Swiz_RF_Temp_RH_data.int = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/switz-sitelocations-17Feb21 copy.csv", na.strings = c("#N/A", "NaN"))

Swiz_RF_Temp_RH_data = Swiz_RF_Temp_RH_data.int[Swiz_RF_Temp_RH_data.int$meanGPP.T<4,]

aov.model = aov(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ species.composition)
aov.model
summary(aov.model)

# Start of Big adjusted all code
# Adjust WUE for Loess graphs
Adjusted_WUE_All = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/ECOSTRESSdata.switzerland.csv")

# match column names
Swiz_RF_Temp_RH_data = Swiz_RF_Temp_RH_data %>% 
  rename(
    Latitude = lat,
    Longitude = long
  )

# merge
str(Adjusted_WUE_All)
str(Swiz_RF_Temp_RH_data)
Adjusted_WUE_All = Adjusted_WUE_All %>% 
  rename(
    site.number = ID 
  )
Big_Adjusted_All = merge(Swiz_RF_Temp_RH_data, Adjusted_WUE_All, by="site.number")
Big_Adjusted_All$M_D = as.character(Big_Adjusted_All$Date)

# create Month_Day Column for plot
Big_Adjusted_All$M_D = substr(Big_Adjusted_All$M_D, 6, 10)

# remove NA from $species.composition and any WUE > 10
Big_Adjusted_All = Big_Adjusted_All[complete.cases(Big_Adjusted_All %>% pull(species.composition)),]
Big_Adjusted_All = Big_Adjusted_All[Big_Adjusted_All$WUE<4,]

# remove all "Beech&Oak&Spruce" columns from $species.composition 
a = which(Big_Adjusted_All$species.composition == "Beech&Oak&Spruce")
Big_Adjusted_All <- Big_Adjusted_All[-a, ]
Big_Adjusted_All$species.composition <- factor(Big_Adjusted_All$species.composition)
Big_Adjusted_All$species.composition = factor(Big_Adjusted_All$species.composition,c("Oak","Beach&Oak","Beech","Beech&Spruce", "Spruce"))

# rename column
Big_Adjusted_All = Big_Adjusted_All %>% 
  rename(
    Species = species.composition
  )

# add space between &, have to change levels for it to work
levels(Big_Adjusted_All$Species) <- c(levels(Big_Adjusted_All$Species), "Beech & Oak")
Big_Adjusted_All$Species[Big_Adjusted_All$Species == 'Beach&Oak'] <- 'Beech & Oak'
levels(Big_Adjusted_All$Species) <- c(levels(Big_Adjusted_All$Species), "Beech & Spruce")
Big_Adjusted_All$Species[Big_Adjusted_All$Species == "Beech&Spruce"] <- "Beech & Spruce"
#Big_Adjusted_All[Big_Adjusted_All$Species %in% "Beech & Spruce",]

# adjust order of levels
Big_Adjusted_All$Species = factor(Big_Adjusted_All$Species, levels = c("Oak", "Beech & Oak", "Beech", "Beech & Spruce", "Spruce"))

# aggregate mean
aggregate(meanGPP.T ~  Species, Big_Adjusted_All, sd)

# read in ET data
ET_June_19 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_June_2019/ET-June-2019-ECO3ETPTJPL-001-results.csv")
ET_June_20 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_June_2020/ET-June-2020-ECO3ETPTJPL-001-results.csv")
ET_July_18 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_July_2018/ET-July-2018-ECO3ETPTJPL-001-results.csv")
ET_July_19 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_July_2019/ET-July-2019-ECO3ETPTJPL-001-results.csv")
ET_July_20 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_July_2020/ET-July-2020-ECO3ETPTJPL-001-results.csv")
ET_Aug_18 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_Aug_2018/ET-Aug-2018-ECO3ETPTJPL-001-results.csv")
ET_Aug_19 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_Aug_2019/ET-Aug-2019-ECO3ETPTJPL-001-results.csv")
ET_Aug_20 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_Aug_2020/ET-Aug-2020-ECO3ETPTJPL-001-results.csv")

# rbind
all_Swiz_ET_data = rbind(ET_June_19, ET_June_20, ET_July_18, ET_July_19, ET_July_20,
                         ET_Aug_18, ET_Aug_19, ET_Aug_20)

colnames(all_Swiz_ET_data)[5] = "Date_total"
colnames(all_Swiz_ET_data)[9] = "ETcanopy"
colnames(all_Swiz_ET_data)[10] = "ETdaily"

# Adjust time zone + separate date and time
all_Swiz_ET_data$date = as_datetime(all_Swiz_ET_data$Date_total)
all_Swiz_ET_data$date.adj = all_Swiz_ET_data$date + hours(1)
all_Swiz_ET_data = separate(data=all_Swiz_ET_data, col=date.adj, into=c("date","time"), sep =" ")

# extract "m-d" from "Y-m-d"
n_last = 5
substr(all_Swiz_ET_data$date, nchar(all_Swiz_ET_data$date) - n_last + 1, nchar(all_Swiz_ET_data$date))
all_Swiz_ET_data$m_d = substr(all_Swiz_ET_data$date, nchar(all_Swiz_ET_data$date) - n_last + 1, nchar(all_Swiz_ET_data$date))

# Make $ID character
all_Swiz_ET_data$ID = as.character(all_Swiz_ET_data$ID)

# seperate into seasonal
Swiz_ET_Data_Summer_18 = all_Swiz_ET_data[all_Swiz_ET_data$date >= "2018-01-01" & all_Swiz_ET_data$date <= "2018-12-31",]
Swiz_ET_Data_Summer_19 = all_Swiz_ET_data[all_Swiz_ET_data$date >= "2019-01-01" & all_Swiz_ET_data$date <= "2019-12-31",]
Swiz_ET_Data_Summer_20 = all_Swiz_ET_data[all_Swiz_ET_data$date >= "2020-01-01" & all_Swiz_ET_data$date <= "2020-12-31",]

# Calculate GPP 
# make coulumns integers to match ET column
all_Swiz_WUE_data$`Orbit Number` = as.integer(all_Swiz_WUE_data$`Orbit Number`)
all_Swiz_WUE_data$`Build ID` = as.integer(all_Swiz_WUE_data$`Build ID`)
all_Swiz_WUE_data$`Scene ID` = as.integer(all_Swiz_WUE_data$`Scene ID`)

# make unique ID by pasting all 8 columns together
all_Swiz_WUE_data$Unique = paste(all_Swiz_WUE_data$Category, all_Swiz_WUE_data$ID, all_Swiz_WUE_data$Latitude,
                                 all_Swiz_WUE_data$Longitude, all_Swiz_WUE_data$Date_total, all_Swiz_WUE_data$`Orbit Number`,
                                 all_Swiz_WUE_data$`Build ID`, all_Swiz_WUE_data$`Scene ID`)

all_Swiz_WUE_data = as.data.frame(all_Swiz_WUE_data)

all_Swiz_ET_data$Unique = paste(all_Swiz_ET_data$Category, all_Swiz_ET_data$ID, all_Swiz_ET_data$Latitude,
                                all_Swiz_ET_data$Longitude, all_Swiz_ET_data$Date_total, all_Swiz_ET_data$Orbit.Number,
                                all_Swiz_ET_data$Build.ID, all_Swiz_ET_data$Scene.ID)

# only need these three columns from ET to make 1 big data frame
Essential_Swiz_ET_data = all_Swiz_ET_data[, c("ETcanopy", "ETdaily", "Unique")]

# merge, final data frame w all data and sites matched up via $unique 
All_Swiz_WUE_ET_Data = merge(all_Swiz_WUE_data, Essential_Swiz_ET_data, by = "Unique")

# WUE (GPP/ET) - no calcs necessary 
head(All_Swiz_WUE_ET_Data$WUEavg)

# WUE (GPP/T) = (WUE x ET Daily) / (ET Daily x ET canopy/100)
All_Swiz_WUE_ET_Data$WUE_GPP_by_T = (All_Swiz_WUE_ET_Data$WUEavg * All_Swiz_WUE_ET_Data$ETdaily) / (All_Swiz_WUE_ET_Data$ETdaily * All_Swiz_WUE_ET_Data$ETcanopy/100)
head(All_Swiz_WUE_ET_Data)

# Remove inf values
All_Swiz_WUE_ET_Data = All_Swiz_WUE_ET_Data %>% 
  filter_if(~is.numeric(.), all_vars(!is.infinite(.)))

# remove values greater than 5
All_Swiz_WUE_ET_Data = All_Swiz_WUE_ET_Data %>% 
  filter(WUE_GPP_by_T <= 5)

# SDS_LST Data
SDS_LST_2018 = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/SDS_LST_Files/Swiz-LST-2018-ECO2LSTE-001-results.csv")
SDS_LST_2019 = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/SDS_LST_Files/Swiz-LST-2019-ECO2LSTE-001-results.csv")
SDS_LST_2020 = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/SDS_LST_Files/Swiz-LST-2020-ECO2LSTE-001-results.csv")

# rbind all data together
all_SDS_LST_data = rbind(SDS_LST_2018, SDS_LST_2019, SDS_LST_2020)

# create $Unique2 to match
All_Swiz_WUE_ET_Data$Unique2 =  paste(All_Swiz_WUE_ET_Data$Latitude, All_Swiz_WUE_ET_Data$Longitude, 
                                      All_Swiz_WUE_ET_Data$Date_total, All_Swiz_WUE_ET_Data$`Orbit Number`,
                                      All_Swiz_WUE_ET_Data$`Build ID`, All_Swiz_WUE_ET_Data$`Scene ID`)

all_SDS_LST_data$Unique2 = paste(all_SDS_LST_data$Latitude, all_SDS_LST_data$Longitude, 
                                 all_SDS_LST_data$Date, all_SDS_LST_data$Orbit.Number,
                                 all_SDS_LST_data$Build.ID, all_SDS_LST_data$Scene.ID)
# get lat long species data
Lat_Long_Species_Data = data.frame(Swiz_RF_Temp_RH_data$Latitude, Swiz_RF_Temp_RH_data$Longitude, 
                                   Swiz_RF_Temp_RH_data$species.composition, Swiz_RF_Temp_RH_data$site.number)
colnames(Lat_Long_Species_Data) = c('Latitude', 'Longitude', 'Species', 'Site Number')

dim(All_Swiz_WUE_ET_Data) #5229
dim(all_SDS_LST_data)     #26336

# merge with WUE/ET, species data
data.new.m2 = merge(All_Swiz_WUE_ET_Data, all_SDS_LST_data, by = 'Unique2') #5168
data.new.m2 = data.new.m2 %>% 
  rename(
    Latitude = Latitude.x,
    Longitude = Longitude.x
  )

# add species from Swiz_RF_Temp_RH_data
data.new.m2$Lat_Long = paste(data.new.m2$Latitude, data.new.m2$Longitude)
Lat_Long_Species_Data$Lat_Long = paste(Lat_Long_Species_Data$Latitude, Lat_Long_Species_Data$Longitude)
data.new = merge(data.new.m2, Lat_Long_Species_Data, by = 'Lat_Long')
data.new$Latitude.y = NULL
data.new$Longitude.y = NULL
data.new$Latitude.y = NULL
data.new$Longitude.y = NULL

# remove NA from $species and any WUE_GPP_by_T > 4
data.new = data.new[complete.cases(data.new %>% pull(Species)),]
data.new = data.new[data.new$WUE_GPP_by_T<4,]

# remove all "Beech&Oak&Spruce" columns from $species.composition 
a = which(data.new$Species == "Beech&Oak&Spruce")
data.new <- data.new[-a, ]
data.new$species.composition <- factor(data.new$Species)
data.new$species.composition = factor(data.new$Species,c("Oak","Beach&Oak","Beech","Beech&Spruce", "Spruce"))
data.new$Species = NULL

# rename column
data.new = data.new %>% 
  rename(
    Species = species.composition
  )

# add space between &, have to change levels for it to work
levels(data.new$Species) <- c(levels(data.new$Species), "Beech & Oak")
data.new$Species[data.new$Species == 'Beach&Oak'] <- 'Beech & Oak'
levels(data.new$Species) <- c(levels(data.new$Species), "Beech & Spruce")
data.new$Species[data.new$Species == "Beech&Spruce"] <- "Beech & Spruce"

# adjust order of levels
data.new$Species = factor(data.new$Species, levels = c("Oak", "Beech & Oak", "Beech", "Beech & Spruce", "Spruce"))

# filter - now say no lowest data
data.WUE.ET.final = data.new[which(data.new$ECO2LSTE_001_SDS_QC_MMD == "0b11" &
                                     data.new$ECO2LSTE_001_SDS_QC_LST_accuracy != c("0b00","0b01") & 
                                     data.new$ECO2LSTE_001_SDS_QC_Mandatory_QA_flags != c("0b11","0b10") & 
                                     data.new$ECO2LSTE_001_SDS_QC_Data_quality_flag != c("0b11","0b10")),]

# remove unnecessary columns/clean columns
data.WUE.ET.final$Unique2 = NULL
data.WUE.ET.final$Unique = NULL
data.WUE.ET.final$Lat_Long = NULL
data.WUE.ET.final = data.WUE.ET.final %>% 
  rename(
    Latitude = Latitude.x,
    Longitude = Longitude.x
  )

# find mean values, SD per site from data.WUE.ET.final, merge new corrected means into data.corrected
test = aggregate(data.WUE.ET.final$WUE_GPP_by_T,list(data.WUE.ET.final$ID),FUN=mean)
test.2 = aggregate(data.WUE.ET.final$WUE_GPP_by_T,list(data.WUE.ET.final$ID),FUN=sd)
names(test) = c("site.number","MeanGPP.T_Corrected")
names(test.2) = c("site.number","SDGPP.T_Corrected")
str(test)
str(data)
data.corrected =  merge(test, test.2) %>%
  merge(data)
str(data.corrected)
str(data.WUE.ET.final)

# what sites got kicked?
kicked.sites = test[which(!test$site.number %in% data.corrected$site.number),]
kicked.sites$site.number

############################## Goldsmith code
setwd("~/Documents/Chapman Research/ECOSTRESS Project/Switzerland")
# need to make sure this data is ok
data.raw = read.csv("ECOSTRESSdata.switzerland-sitelocations-17Feb21.csv",header=TRUE,na.strings="#N/A")
data = data.raw[which(data.raw$meanGPP.T < 4),]

# add space around "&"
levels(data$species.composition) <- c(levels(data$species.composition), "Beech & Oak")
data$species.composition[data$species.composition == 'Beech&Oak'] <- 'Beech & Oak'
levels(data$species.composition) <- c(levels(data$species.composition), "Beech & Spruce")
data$species.composition[data$species.composition == "Beech&Spruce"] <- "Beech & Spruce"

# remove excess factors and NA
data$species.composition <- factor(data$species.composition)
data$species.composition = factor(data$species.composition,c("Oak","Beech & Oak","Beech","Beech & Spruce", "Spruce"))
data = data[complete.cases(data %>% pull(species.composition)),]

#rename
data = data %>% 
  rename(
    Species = species.composition
  )
# remove the extra site 1083
data = data[-c(62),]
dim(data)
length(unique(data$site.number))

#### Look up how to test linear model fits in R: 
# Pearson
cor.test(data$meanGPP.T, data$mean.annual.P)

# aggregate mean
aggregate(meanGPP.T ~  Species, data, mean)

Oak.data = data[data$species.composition == "Oak",]
Oak.data.b = na.omit(Oak.data$mean.annual.RH)
Beech.data = data[data$species.composition == "Beech",]
Beech.data.b = na.omit(Beech.data$mean.annual.RH)
spruce.data = data[data$species.composition == "Spruce",]
spruce.data.b = na.omit(spruce.data$mean.annual.RH)
mean(spruce.data.b)

# tukey posthoc test
wue.lm <- lm(meanGPP.T ~ species.composition, data = data)
wue.av <- aov(wue.lm)
summary(wue.av)
tukey.test <- TukeyHSD(wue.av)
tukey.test

################################ Nutrients
# Nutrients Data Set
install.packages("readxl")
library("readxl")
library(tidyverse)
nutrient.data = as.data.frame(read_excel("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/nutrientsbufiei1519 (1).xls"))
nutrient.data.2019 = nutrient.data[nutrient.data$JAHR == 2019,]
nutrient.data.2019$STNRNEU = as.integer(nutrient.data.2019$STNRNEU)
nutrient.data.2019$`SPECIES$` = as.factor(nutrient.data.2019$`SPECIES$`)
nutrient.data.2019 = nutrient.data.2019 %>% 
  rename(
    site.number = STNRNEU,
    species.composition = `SPECIES$`
  )
data = data[-c(76),]

# merge
nutrient.data.2019.agg = aggregate(nutrient.data.2019$STICKST, list(nutrient.data.2019$site.number), FUN = mean, na.rm = TRUE)
names(nutrient.data.2019.agg) = c("site.number", "Leaf.Nitrogen")
nutrient.data.2019.agg$Calcium = aggregate(nutrient.data.2019$CALCIUM, list(nutrient.data.2019$site.number), FUN = mean, na.rm = TRUE)

all.data = merge(data, nutrient.data.2019.agg, by='site.number')
str(all.data)

#try w/ 3 year data first
# when is N data from, same graphs for Ozone try with just 2019 WUE if possible
#for all enviro - ozone, n dep and leaf nutrients
#for all color by species comp
# plot N dep x axis by leaf N y axis.

#data.WUE.ET.final = individual site data
#data.corrected = filtered, mean data
unique(data.corrected$Species)

# read in tables as data frames
nutrientsbufiei1519 = as.data.frame(read_xls("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Environmental_Factor_Data/nutrientsbufiei1519 (1).xls"))
ozon1519 = as.data.frame(read_xls("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Environmental_Factor_Data/ozon1519.xls"))

# rename columns
nutrientsbufiei1519 = nutrientsbufiei1519 %>% 
  rename(
    year = JAHR, site.number = STNRNEU, nitrogen = STICKST, phosphorus = PHOSPHOR,
    potassium = KALIUM, calcium = CALCIUM, magnesium = MAGNES, manganese = MANGAN,
    NP.ratio = NP, NK.ratio = NK, NMg.ratio = NMG, Species = "SPECIES$"
  )

#nutrientsbufiei1519.col.names = colnames(nutrientsbufiei1519[,3:11])
#test.3 = c(3:11)
#means_list = list()
#for (i in test.3) {
 # new.data = nutrientsbufiei1519[c(2,i)]
  #colnames(new.data) = c("id", "value")
  #Ag <- aggregate(value ~ id, mean,data=new.data)
  #means_list[[i]] = Ag
  #id = new.data[c(2)]
  #means_list[[i]] = aggregate(nutrientsbufiei1519$i,list(nutrientsbufiei1519$site.number),FUN=mean)
#}
# could create 3rd column, that IDs what each value is from ex. nitroen, phos, etc...

# take means and rename columns, then merge
nutrientsbufiei1519.nitrogen = aggregate(nutrientsbufiei1519$nitrogen,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.nitrogen) = c("site.number", "mean.nitrogen")
nutrientsbufiei1519.phosphorus = aggregate(nutrientsbufiei1519$phosphorus,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.phosphorus) = c("site.number", "mean.phosphorus")
nutrientsbufiei1519.potassium = aggregate(nutrientsbufiei1519$potassium,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.potassium) = c("site.number", "mean.potassium")
nutrientsbufiei1519.calcium = aggregate(nutrientsbufiei1519$calcium,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.calcium) = c("site.number", "mean.calcium")
nutrientsbufiei1519.magnesium = aggregate(nutrientsbufiei1519$magnesium,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.magnesium) = c("site.number", "mean.magnesium")
nutrientsbufiei1519.manganese = aggregate(nutrientsbufiei1519$manganese,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.manganese) = c("site.number", "mean.manganese")
nutrientsbufiei1519.NP.ratio = aggregate(nutrientsbufiei1519$NP.ratio,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.NP.ratio) = c("site.number", "mean.nNP.ratio")
nutrientsbufiei1519.NK.ratio = aggregate(nutrientsbufiei1519$NK.ratio,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.NK.ratio) = c("site.number", "mean.NK.ratio")
nutrientsbufiei1519.NMg.ratio = aggregate(nutrientsbufiei1519$NMg.ratio,list(nutrientsbufiei1519$site.number),FUN=mean)
names(nutrientsbufiei1519.NMg.ratio) = c("site.number", "mean.NMg.ratio")

nutrientsbufiei1519_means = merge(nutrientsbufiei1519.nitrogen, nutrientsbufiei1519.phosphorus) %>%
  merge(nutrientsbufiei1519.potassium) %>%
  merge(nutrientsbufiei1519.calcium) %>%
  merge(nutrientsbufiei1519.magnesium) %>%
  merge(nutrientsbufiei1519.manganese) %>%
  merge(nutrientsbufiei1519.NP.ratio) %>%
  merge(nutrientsbufiei1519.NK.ratio) %>%
  merge(nutrientsbufiei1519.NMg.ratio)

# take 2019 data only
WUE_ET_Data_2019 = All_Swiz_WUE_ET_Data[All_Swiz_WUE_ET_Data$date >= "2019-01-01" & All_Swiz_WUE_ET_Data$date <= "2019-12-31",]
WUE_ET_Data_2019 = WUE_ET_Data_2019 %>% 
  rename(
    site.number = ID
  )


# subset WUE to only means of WUE per site for 2019 data
WUE_Data_2019_Means = aggregate(WUE_ET_Data_2019$WUE_GPP_by_T,list(WUE_ET_Data_2019$site.number),FUN=mean)
names(WUE_Data_2019_Means) = c("site.number", "meanGPP.T")
site_species = data.frame(data.WUE.ET.final$`Site Number`, data.WUE.ET.final$Species)
site_species = site_species %>% 
  rename(
    site.number = data.WUE.ET.final..Site.Number.,
    Species = data.WUE.ET.final.Species
  )

data.WUE.2019.site.means = unique(merge(WUE_Data_2019_Means, site_species))

# 2019 nutrients data, ceate data frame then need to merge this with data.WUE.2019.site.means
nutrientsbufiei1519_2019 = nutrientsbufiei1519[which(nutrientsbufiei1519$year == "2019"),]
nutrient.2019.sitemeans_calcium = aggregate(nutrientsbufiei1519_2019$calcium,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_calcium) = c("site.number", "mean.calcium")
nutrient.2019.sitemeans_nitrogen = aggregate(nutrientsbufiei1519_2019$nitrogen,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_nitrogen) = c("site.number", "mean.nitrogen")
nutrient.2019.sitemeans_phosphorus = aggregate(nutrientsbufiei1519_2019$phosphorus,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_phosphorus) = c("site.number", "mean.phosphorus")
nutrient.2019.sitemeans_potassium = aggregate(nutrientsbufiei1519_2019$potassium,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_potassium) = c("site.number", "mean.potassium")
nutrient.2019.sitemeans_magnesium = aggregate(nutrientsbufiei1519_2019$magnesium,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_magnesium) = c("site.number", "mean.magnesium")
nutrient.2019.sitemeans_manganese = aggregate(nutrientsbufiei1519_2019$manganese,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_manganese) = c("site.number", "mean.manganese")
nutrient.2019.sitemeans_NP.ratio = aggregate(nutrientsbufiei1519_2019$NP.ratio,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_NP.ratio) = c("site.number", "mean.NP.ratio")
nutrient.2019.sitemeans_NK.ratio = aggregate(nutrientsbufiei1519_2019$NK.ratio,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_NK.ratio) = c("site.number", "mean.NK.ratio")
nutrient.2019.sitemeans_NMg.ratio = aggregate(nutrientsbufiei1519_2019$NMg.ratio,list(nutrientsbufiei1519_2019$site.number),FUN=mean,na.rm=T)
names(nutrient.2019.sitemeans_NMg.ratio) = c("site.number", "mean.NMg.ratio")

nutrient.2019.sitemeans = merge(nutrient.2019.sitemeans_calcium, nutrient.2019.sitemeans_nitrogen) %>%
  merge(nutrient.2019.sitemeans_phosphorus) %>%
  merge(nutrient.2019.sitemeans_potassium) %>%
  merge(nutrient.2019.sitemeans_magnesium) %>%
  merge(nutrient.2019.sitemeans_manganese) %>%
  merge(nutrient.2019.sitemeans_NP.ratio) %>%
  merge(nutrient.2019.sitemeans_NK.ratio) %>%
  merge(nutrient.2019.sitemeans_NMg.ratio)

# merge data.WUE.2019.site.means with nutrient.2019.sitemeans
WUE.2019.Nutrient.2019.data = merge(data.WUE.2019.site.means, nutrient.2019.sitemeans, by = "site.number", na.rm=T)

# now for all WUE means, not just 2019 
# rename for all data
All_Swiz_WUE_ET_Data = All_Swiz_WUE_ET_Data %>% 
  rename(
    site.number = ID
  )

# subset WUE to  means of WUE per site for all data
WUE_Data_Means = aggregate(All_Swiz_WUE_ET_Data$WUE_GPP_by_T,list(All_Swiz_WUE_ET_Data$site.number),FUN=mean)
names(WUE_Data_Means) = c("site.number", "meanGPP.T")
site_species = data.frame(data.WUE.ET.final$`Site Number`, data.WUE.ET.final$Species)
site_species = site_species %>% 
  rename(
    site.number = data.WUE.ET.final..Site.Number.,
    Species = data.WUE.ET.final.Species
  )

data.all.WUE.site.means = unique(merge(WUE_Data_Means, site_species))

# merge data.all.WUE.site.means with nutrient.2019.sitemeans
WUE.Nutrient.2019.data = merge(data.all.WUE.site.means, nutrient.2019.sitemeans, by = "site.number", na.rm=T)




# ozone
ozon1519 = ozon1519 %>% 
  rename(
    site.number = STNRNEU, year = JAHR, ozone.flux.beech = BD1SW, ozone.flux.spruce = FD1SW
  )

# get means for each site at each year
ozon1519$site.avg <- rowMeans(ozon1519[,c('ozone.flux.beech', 'ozone.flux.spruce')], na.rm=TRUE)

#ozone1519.beech = aggregate(ozon1519$ozone.flux.beech,list(ozon1519$site.number),FUN=mean)
#ozone1519.spruce = aggregate(ozon1519$ozone.flux.spruce,list(ozon1519$site.number),FUN=mean)
#names(ozone1519.beech) = c("site.number", "mean.ozone.flux.beech")
#names(ozone1519.spruce) = c("site.number", "mean.ozone.flux.spruce")
#ozone1519_means = merge(ozone1519.beech,ozone1519.spruce)

# take 2019 Ozone data only, merge with other 2019 data
ozon1519_2019 = ozon1519[which(ozon1519$year == "2019"),]
ozon1519_2019_2 = unique(merge(ozon1519_2019, site_species))
GPP.T_2019_Ozone_2019_data = merge(ozon1519_2019_2, WUE_Data_2019_Means, by = "site.number")

# mean of all ozone data merge with mean WUE data
ozone.all.sitemeans = aggregate(ozon1519$site.avg,list(ozon1519$site.number),FUN=mean,na.rm=T)
names(ozone.all.sitemeans) = c("site.number", "mean.site.avg")
ozon1519_2 = unique(merge(ozone.all.sitemeans, site_species))
GPP.T_Ozone_data = merge(ozon1519_2, WUE_Data_Means, by = "site.number")

# Read in carbon 13 validation data
carbon13validationdata = as.data.frame(read_csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Environmental_Factor_Data/carbon13validationdata.csv"))
carbon13validationdata$X1 = NULL
carbon13validationdata = carbon13validationdata %>% 
  rename(
    site.number = siteID
    )
# merge with site mean WUE
carbon13validationdata_wue = merge(carbon13validationdata, WUE_Data_2019_Means, by = "site.number")
carbon13validationdata_final = unique(merge(carbon13validationdata_wue, site_species))

# ADDITIONAL POINTS
New_Point_1 = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Extra_Sites/Site_1/WUE/New-Point-18-19-20-ECO4WUE-001-results (1).csv")
New_Point_1_SDS_LST = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Extra_Sites/Site_1/SDS_LST/New-Point-SDS-LST-ECO2LSTE-001-results.csv")
New_Point_2 = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Extra_Sites/Site_2/WUE/Second-New-Point-18-19-20-ECO4WUE-001-results.csv")
New_Point_2_SDS_LST = read.csv("/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Extra_Sites/Site_2/SDS_LST/Second-Point-SDS-LST-ECO2LSTE-001-results.csv")

# label each site
New_Point_1$Site = "Site_1"
New_Point_2$Site = "Site_2"
# bind
All_New_Point_WUE_Data = rbind(New_Point_1, New_Point_2)
All_New_Point_LST_SDS_Data = rbind(New_Point_1_SDS_LST, New_Point_2_SDS_LST)

# create $Unique
All_New_Point_WUE_Data$Unique = paste(All_New_Point_WUE_Data$Latitude, All_New_Point_WUE_Data$Longitude, 
                                      All_New_Point_WUE_Data$Date, All_New_Point_WUE_Data$Orbit.Number,
                                      All_New_Point_WUE_Data$Scene.ID, All_New_Point_WUE_Data$Build.ID)

All_New_Point_LST_SDS_Data$Unique = paste(All_New_Point_LST_SDS_Data$Latitude, All_New_Point_LST_SDS_Data$Longitude, 
                                           All_New_Point_LST_SDS_Data$Date, All_New_Point_LST_SDS_Data$Orbit.Number,
                                           All_New_Point_LST_SDS_Data$Scene.ID, All_New_Point_LST_SDS_Data$Build.ID)

# merge
All_New_Point_Data = merge(All_New_Point_WUE_Data, All_New_Point_LST_SDS_Data, by = "Unique")

# filter by SDS_LST
All_New_Point_Data_filtered = All_New_Point_Data[which(All_New_Point_Data$ECO2LSTE_001_SDS_QC_MMD == "0b11" &
                                                 All_New_Point_Data$ECO2LSTE_001_SDS_QC_LST_accuracy != c("0b00","0b01") & 
                                                 All_New_Point_Data$ECO2LSTE_001_SDS_QC_Mandatory_QA_flags != c("0b11","0b10") & 
                                                 All_New_Point_Data$ECO2LSTE_001_SDS_QC_Data_quality_flag != c("0b11","0b10")),]

# Filter by WUE < 4
All_New_Point_Data_final = All_New_Point_Data_filtered[which(All_New_Point_Data_filtered$ECO4WUE_001_Water_Use_Efficiency_WUEavg <= 4),]

All_New_Point_Data_final$ECO4WUE_001_Water_Use_Efficiency_WUEavg
Site_1 = All_New_Point_Data_filtered[which(All_New_Point_Data_filtered$Site == "Site_1"),]
Site_2 = All_New_Point_Data_filtered[which(All_New_Point_Data_filtered$Site == "Site_2"),]
sd(Site_2$ECO4WUE_001_Water_Use_Efficiency_WUEavg)

######################################################################################################
######################################################################################################
######################################################################################################

# PLOTS

# env factor plots, not used this way anymore
par(mfrow = c(3,1))
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ mean.annual.P)
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ mean.annual.T)
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ mean.annual.RH)

par(mfrow = c(3,1))
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ elevation)
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ slope.degree)
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ aspect.degree)

par(mfrow = c(3,1))
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ total.n.deposition)
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ oxidized.nitrogen)
plot(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ reduced.nitrogen)
abline(lm(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ reduced.nitrogen))
summary(lm(data = Swiz_RF_Temp_RH_data, meanGPP.T ~ reduced.nitrogen))

# boxplot big adjusted
boxplot(data = Big_Adjusted_All, meanGPP.T ~ Species, frame.plot = FALSE,
        ylab = "WUE (GPP/T)", xlab = "Dominant Tree Species", cex.lab=1.75, cex.axis=1.75)
box(bty="l")

# some LOESS
# colored by species
str(data.WUE.ET.final)

# rename for correct legend
data.WUE.ET.final = data.WUE.ET.final %>% 
  rename(
    'Species Composition' = Species,
  )

Swiz_loess_Adjusted_WUE_All = ggplot(data.WUE.ET.final, aes(x = as.POSIXct(m_d,format = "%m-%d"), 
                                                            y = WUE_GPP_by_T, color = `Species Composition`)) +
  geom_point(alpha = 0.5) + geom_smooth(se = FALSE) + theme_bw() + ggtitle(" ") +
  theme(plot.title = element_text(size = 12, family = "Tahoma"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))  +
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)"))) 

# All Summers 
Swiz_loess_WUE_All = ggplot(data.WUE.ET.final, aes(x = as.POSIXct(m_d,format = "%m-%d"), 
                                                   y = WUEavg, color = ID)) + 
  geom_point(alpha = 0.5) + geom_smooth(se = FALSE) + theme_bw() + ggtitle(" ") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))+
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2018 WUE
Swiz_loess_WUE_Summer_18 = ggplot(Swiz_WUE_Data_Summer_18, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2018") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2019 WUE
Swiz_loess_WUE_Summer_19 = ggplot(Swiz_WUE_Data_Summer_19, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2019") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2020 WUE
Swiz_loess_WUE_Summer_20 = ggplot(Swiz_WUE_Data_Summer_20, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2020") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# ggplot geom_density histogram of WUE (GPP/ET) compared to WUE (GPP/T)
# WUE (GPP/T)
Plot_Swiz_WUE_GPP_T = ggplot(data.WUE.ET.final, aes(x=WUE_GPP_by_T)) + 
  geom_histogram(aes(y=..density..), binwidth=.05,color="#87CEFA", fill="white") +
  geom_density(alpha=.2, fill="#87CEFA") + xlim(0,4) + ylim(0,1.25) + 
  geom_vline(aes(xintercept=mean(WUE_GPP_by_T)), color="blue", linetype="dashed", size=1) +
  xlab(expression(paste("WUE (GPP/T) (g C kg"^"-1 ",H[2],"O)"))) + ylab("Density") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# WUE (GPP/ET)
Plot_Swiz_WUE_GPP_ET = ggplot(data.WUE.ET.final, aes(x=WUEavg)) + 
  geom_histogram(aes(y=..density..), binwidth=.05,color="#F08080", fill="white") +
  geom_density(alpha=.2, fill="#F08080") + xlim(0,4) + ylim(0,1.25) + 
  geom_vline(aes(xintercept=mean(WUEavg)), color="red", linetype="dashed", size=1) +
  xlab(expression(paste("WUE (GPP/ET) (g C kg"^"-1 ",H[2],"O)"))) + ylab("Density") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

mean(data.WUE.ET.final$WUE_GPP_by_T)

# 2 Panel Image
Plot_Swiz_Both_WUE = ggarrange(Plot_Swiz_WUE_GPP_ET, Plot_Swiz_WUE_GPP_T,
                               labels = c("A)", "B)"),
                               font.label = list(size = 12, color = "black"),
                               ncol = 1, nrow = 2)

##### 9 panel using ggplot

# elevation
GPP.T_Elevation = data.corrected %>%
  ggplot(aes(elevation,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x="Elevation (m asl)")+ theme_bw() +
  theme(plot.title = element_text(size = 11, family = "Tahoma"),
        text = element_text(size = 11, family = "Tahoma"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11), legend.position="none") + 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95)

# lm use * first, if not sig then use +...
install.packages("lme4")
library("lme4")
model_WUE_Elevation = lm(MeanGPP.T_Corrected ~ Species + elevation, data = data.corrected)
model_WUE_Elevation_lmer = lmer(MeanGPP.T_Corrected ~ elevation + (1|Species), data = data.corrected)
summary(model_WUE_Elevation)
aov(model_WUE_Elevation)

# slope
GPP.T_Slope = data.corrected %>%
  ggplot(aes(slope.degree,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x="Slope (\u00B0)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95)

# lm
model_WUE_Slope = lm(MeanGPP.T_Corrected ~ Species + slope.degree, data = data.corrected)
model_WUE_Slope_lmer = lmer(MeanGPP.T_Corrected ~ slope.degree + (1|Species), data = data.corrected)
summary(model_WUE_Slope)
aov(model_WUE_Slope)

# aspect
GPP.T_Aspect = data.corrected %>%
  ggplot(aes(aspect.degree,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x="Aspect (\u00B0)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95)

# lm
model_WUE_Aspect = lm(MeanGPP.T_Corrected ~ Species + aspect.degree, data = data.corrected)
model_WUE_Aspect_lmer = lmer(MeanGPP.T_Corrected ~ aspect.degree + (1|Species), data = data.corrected)
summary(model_WUE_Aspect)
aov(model_WUE_Aspect)

# temp
GPP.T_Temp = data.corrected %>%
  ggplot(aes(mean.annual.T,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x="Mean annual temperature (\u00B0C)")+ theme_bw() +
  theme(plot.title = element_text(size = 11, family = "Tahoma"),
        text = element_text(size = 11, family = "Tahoma"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11), legend.position="none") + 
  stat_cor(method = "pearson", label.y.npc = 0.95)

# lm
model_WUE_Temp = lm(MeanGPP.T_Corrected ~ Species * mean.annual.T, data = data.corrected)
model_WUE_Temp_lmer = lmer(MeanGPP.T_Corrected ~ mean.annual.T + (1|Species), data = data.corrected)
summary(model_WUE_Temp)
aov(model_WUE_Temp)

# precip
GPP.T_Precipitation = data.corrected %>%
  ggplot(aes(mean.annual.P,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x="Mean annual precipitation (mm)")+ theme_bw() +
  theme(plot.title = element_text(size = 11, family = "Tahoma"),
        text = element_text(size = 11, family = "Tahoma"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11), legend.position="none") + 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95)

# lm
model_WUE_Precip = lm(MeanGPP.T_Corrected ~ Species + mean.annual.P, data = data.corrected)
model_WUE_Precip_lmer = lmer(MeanGPP.T_Corrected ~ mean.annual.P + (1|Species), data = data.corrected)
summary(model_WUE_Precip)
aov(model_WUE_Precip)

# Rel. Humidity
GPP.T_Rel_Humidity = data.corrected %>%
  ggplot(aes(mean.annual.RH,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x="Mean annual relative humidity (%)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.y.npc = 0.95)

# lm
model_WUE_RH = lm(MeanGPP.T_Corrected ~ Species + mean.annual.RH, data = data.corrected)
model_WUE_RH_lmer = lmer(MeanGPP.T_Corrected ~ mean.annual.RH + (1|Species), data = data.corrected)
summary(model_WUE_RH)
aov(model_WUE_RH)

# oxidized n
GPP.T_Oxidized_N = data.corrected %>%
  ggplot(aes(oxidized.nitrogen,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x=expression(paste(NO[3],""^"-"," (kg N ",ha^-1, yr^-1,")")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95)

# lm
model_WUE_Oxidized_N = lm(MeanGPP.T_Corrected ~ Species + oxidized.nitrogen, data = data.corrected)
model_WUE_Oxidized_N_lmer = lmer(MeanGPP.T_Corrected ~ oxidized.nitrogen + (1|Species), data = data.corrected)
summary(model_WUE_Oxidized_N)
aov(model_WUE_Oxidized_N)

# reduced n
GPP.T_Reduced_N = data.corrected %>%
  ggplot(aes(reduced.nitrogen,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x=expression(paste(NH[4],""^"+"," (kg N ",ha^-1, yr^-1,")")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Reduced_N = lm(MeanGPP.T_Corrected ~ Species * reduced.nitrogen, data = data.corrected)
model_WUE_Reduced_N_lmer = lmer(MeanGPP.T_Corrected ~ reduced.nitrogen + (1|Species), data = data.corrected)
summary(model_WUE_Reduced_N)
aov(model_WUE_Reduced_N)

# n dep
GPP.T_N_Deposition = data.corrected %>%
  ggplot(aes(total.n.deposition,MeanGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), x=expression(paste("Total nitrogen deposition (kg N ",ha^-1, yr^-1,")")))+ theme_bw() +
  theme(plot.title = element_text(size = 11, family = "Tahoma"),
        text = element_text(size = 11, family = "Tahoma"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_N_Deposition = lm(MeanGPP.T_Corrected ~ Species + total.n.deposition, data = data.corrected)
model_WUE_N_Deposition_lmer = lmer(MeanGPP.T_Corrected ~ total.n.deposition + (1|Species), data = data.corrected)
summary(model_WUE_N_Deposition)
aov(model_WUE_N_Deposition)

# Panel images
EF_WUE_6Panel = ggarrange(GPP.T_Elevation, GPP.T_Slope, GPP.T_Aspect,
                          GPP.T_Temp, GPP.T_Precipitation, GPP.T_Rel_Humidity,
                          labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
                          font.label = list(size = 12),
                          ncol = 3, nrow = 2)

Nitrogen_WUE_4Panel = ggarrange(GPP.T_Oxidized_N, GPP.T_Reduced_N, GPP.T_N_Deposition,
                                GPP.T_Leaf_Nitrogen, 
                                labels = c("A)", "B)", "C)", "D)"),
                                font.label = list(size = 12),
                                ncol = 2, nrow = 2)

###### 9 panel SD

# elevation
GPP.T_Elevation_SD = data.corrected %>%
  ggplot(aes(elevation,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x="Elevation (m asl)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Elevation_SD = lm(SDGPP.T_Corrected ~ Species + elevation, data = data.corrected)
model_WUE_Elevation_SD_lmer = lmer(SDGPP.T_Corrected ~ elevation + (1|Species), data = data.corrected)
summary(model_WUE_Elevation)
aov(model_WUE_Elevation)

# slope
GPP.T_Slope_SD = data.corrected %>%
  ggplot(aes(slope.degree,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x="Slope (\u00B0)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Slope_SD = lm(SDGPP.T_Corrected ~ Species + slope.degree, data = data.corrected)
model_WUE_Slope_SD_lmer = lmer(SDGPP.T_Corrected ~ slope.degree + (1|Species), data = data.corrected)
summary(model_WUE_Slope_SD)
aov(model_WUE_Slope_SD)

# aspect
GPP.T_Aspect_SD = data.corrected %>%
  ggplot(aes(aspect.degree,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x="Aspect (\u00B0)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Aspect_SD = lm(SDGPP.T_Corrected ~ Species + aspect.degree, data = data.corrected)
model_WUE_Aspect_SD_lmer = lmer(SDGPP.T_Corrected ~ aspect.degree + (1|Species), data = data.corrected)
summary(model_WUE_Aspect_SD)
aov(model_WUE_Aspect_SD)

# temp
GPP.T_Temp_SD = data.corrected %>%
  ggplot(aes(mean.annual.T,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x="Mean annual temperature (\u00B0C)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Temp_SD = lm(SDGPP.T_Corrected ~ Species + mean.annual.T, data = data.corrected)
model_WUE_Temp_SD_lmer = lmer(SDGPP.T_Corrected ~ mean.annual.T + (1|Species), data = data.corrected)
summary(model_WUE_Temp_SD)
aov(model_WUE_Temp_SD)

# precip
GPP.T_Precipitation_SD = data.corrected %>%
  ggplot(aes(mean.annual.P,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x="Mean annual precipitation (mm)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Precipitation_SD = lm(SDGPP.T_Corrected ~ Species + mean.annual.P, data = data.corrected)
model_WUE_Precipitation_SD_lmer = lmer(SDGPP.T_Corrected ~ mean.annual.P + (1|Species), data = data.corrected)
summary(model_WUE_Precipitation_SD)
aov(model_WUE_Precipitation_SD)

# Rel. Humidity
GPP.T_Rel_Humidity_SD = data.corrected %>%
  ggplot(aes(mean.annual.RH,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x="Mean annual relative humidity (%)")+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Rel_Humidity_SD = lm(SDGPP.T_Corrected ~ Species + mean.annual.RH, data = data.corrected)
model_WUE_Rel_Humidity_SD_lmer = lmer(SDGPP.T_Corrected ~ mean.annual.RH + (1|Species), data = data.corrected)
summary(model_WUE_Rel_Humidity_SD)
aov(model_WUE_Rel_Humidity_SD)

# oxidized n
GPP.T_Oxidized_N_SD = data.corrected %>%
  ggplot(aes(oxidized.nitrogen,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x=expression(paste(NO[3],""^"-"," (kg N ",ha^-1, yr^-1,")")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Oxidized_N_SD = lm(SDGPP.T_Corrected ~ Species + oxidized.nitrogen, data = data.corrected)
model_WUE_Oxidized_N_SD_lmer = lmer(SDGPP.T_Corrected ~ oxidized.nitrogen + (1|Species), data = data.corrected)
summary(model_WUE_Oxidized_N_SD)
aov(model_WUE_Oxidized_N_SD)

# reduced n
GPP.T_Reduced_N_SD = data.corrected %>%
  ggplot(aes(reduced.nitrogen,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x=expression(paste(NH[4],""^"+"," (kg N ",ha^-1, yr^-1,")")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_Reduced_N_SD = lm(SDGPP.T_Corrected ~ Species + reduced.nitrogen, data = data.corrected)
model_WUE_Reduced_N_SD_lmer = lmer(SDGPP.T_Corrected ~ reduced.nitrogen + (1|Species), data = data.corrected)
summary(model_WUE_Reduced_N_SD)
aov(model_WUE_Reduced_N_SD)

# n dep
GPP.T_N_Deposition_SD = data.corrected %>%
  ggplot(aes(total.n.deposition,SDGPP.T_Corrected)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("SD WUE (g C ",kg^-1," ",H[2],O,")")), x=expression(paste("Total nitrogen deposition (kg N ",ha^-1, yr^-1,")")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm
model_WUE_N_Deposition_SD = lm(SDGPP.T_Corrected ~ Species + total.n.deposition, data = data.corrected)
model_WUE_N_Deposition_SD_lmer = lmer(SDGPP.T_Corrected ~ total.n.deposition + (1|Species), data = data.corrected)
summary(model_WUE_N_Deposition_SD)
aov(model_WUE_N_Deposition_SD)

Plot_Swiz_Both_WUE_SD = ggarrange(GPP.T_Elevation_SD, GPP.T_Slope_SD, GPP.T_Aspect_SD, GPP.T_Temp_SD, 
                                  GPP.T_Precipitation_SD, GPP.T_Rel_Humidity_SD, GPP.T_Oxidized_N_SD, 
                                  GPP.T_Reduced_N_SD, GPP.T_N_Deposition_SD,
                                  ncol = 3, nrow = 3)

###### wue boxplot by species

wue_boxplot = ggplot(data.WUE.ET.final, aes(x = `Species Composition`, y = WUE_GPP_by_T, fill = `Species Composition`)) + 
  geom_boxplot(notch=TRUE) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Species Composition")))+ theme_bw() +
  theme(plot.title = element_text(size = 12, family = "Tahoma"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), legend.position="none")

ggsave(wue_boxplot, plot = wue_boxplot, device = "pdf",
       path = "/Users/brandonbernardo/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Final_Plots"
      )

beech_spruce_ex = data.WUE.ET.final[which(data.WUE.ET.final$`Species Composition` == "Beech & Spruce"),]
sd(spruce_ex$WUE_GPP_by_T)

# boxplot stats
ggplot_build(wue_boxplot)$data

# ANOVA test
res_aov = aov(WUE_GPP_by_T ~ `Species Composition`, data = data.WUE.ET.final)
summary(res_aov)

# histogram
hist(res_aov$residuals)

# QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# 2 panel image loess + boxplot
WUE_LOESS_Boxplot = ggarrange(Swiz_loess_Adjusted_WUE_All, wue_boxplot,
          labels = c("A)", "B)"),
          font.label = list(size = 12, color = "black"),
          ncol = 2, nrow = 1)

# Leaf Nitrogen
GPP.T_Leaf_Nitrogen = all.data %>%
  ggplot(aes(Leaf.Nitrogen,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Leaf Nitrogen")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.text=element_text(size=15),
        legend.title=element_text(size=15), legend.position="none") + 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.accuracy = 0.001, r.accuracy = 0.001)

# lm
model_WUE_Leaf_Nitrogen = lm(meanGPP.T ~ Species * Leaf.Nitrogen, data = all.data)
model_WUE_Leaf_Nitrogen_lmer = lmer(meanGPP.T ~ Leaf.Nitrogen + (1|Species), data = all.data)
summary(model_WUE_Leaf_Nitrogen)
aov(model_WUE_Leaf_Nitrogen)

# Calcium
GPP.T_Calcium = all.data %>%
  ggplot(aes(Calcium,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Calcium")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")
all.data$Calcium

# lm - need?
model_WUE_Leaf_Calcium = lm(meanGPP.T ~ Species * Calcium, data = all.data)
model_WUE_Leaf_Calcium_lmer = lmer(meanGPP.T ~ Calcium + (1|Species), data = all.data)
summary(model_WUE_Leaf_Calcium)
aov(model_WUE_Leaf_Calcium)

# Requested Plots
# 1
# Mean leaf N x N dep
N_Dep_Mean_Nitrogen = data.corrected %>%
  ggplot(aes(total.n.deposition,mean.nitrogen)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean Nitrogen")), 
       x=expression(paste("Nitrogen Deposition"))) + theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# lm - Need?
model_N_Leaf_Calcium = lm(mean.nitrogen ~ Species * total.n.deposition, data = data.corrected)
model_N_Leaf_Calcium_lmer = lmer(mean.nitrogen ~ total.n.deposition + (1|Species), data = data.corrected)
summary(model_N_Leaf_Calcium)
aov(model_N_Leaf_Calcium)

# Mean leaf N x ammonium - MISSING
N_Dep_Mean_Nitrogen = data.corrected %>%
  ggplot(aes(total.n.deposition,mean.nitrogen)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean Nitrogen")), 
       x=expression(paste("Nitrogen Deposition"))) + theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# Mean leaf N x Nitrate - MISSING
N_Dep_Mean_Nitrogen = data.corrected %>%
  ggplot(aes(total.n.deposition,mean.nitrogen)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean Nitrogen")), 
       x=expression(paste("Nitrogen Deposition"))) + theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# 2
# Mean WUE x mean ozone - not by species
GPP.T_Ozone = GPP.T_Ozone_data %>%
  ggplot(aes(mean.site.avg,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Ozone Site Averages")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# lm 
model_WUE_Leaf_Ozone = lm(meanGPP.T ~ Species + mean.site.avg, data = GPP.T_Ozone_data)
model_WUE_Leaf_Ozone_lmer = lmer(meanGPP.T ~ mean.site.avg + (1|Species), data = GPP.T_Ozone_data)
summary(model_WUE_Leaf_Ozone)
aov(model_WUE_Leaf_Ozone)

# 3 
# 2019 WUE x 2019 mean ozone - not by species
GPP.T_2019_Ozone_2019 = GPP.T_2019_Ozone_2019_data %>%
  ggplot(aes(site.avg,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean 2019 WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("2019 Ozone Site Averages")))+ theme_bw() +
  theme(plot.title = element_text(size = 15, family = "Tahoma"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15), legend.position="none")+ 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.digits = signif(2), r.digits = signif(2)) 

# lm 
model_WUE2019_Leaf_Ozone2019 = lm(meanGPP.T ~ Species + site.avg, data = GPP.T_2019_Ozone_2019_data)
model_WUE2019_Leaf_Ozone2019_lmer = lmer(meanGPP.T ~ site.avg + (1|Species), data = GPP.T_2019_Ozone_2019_data)
summary(model_WUE2019_Leaf_Ozone2019)
aov(model_WUE2019_Leaf_Ozone2019)

# 4 
# 2019 WUE x 2019 Leaf N
GPP.T_2019_Nitrogen_2019 = WUE.2019.Nutrient.2019.data %>%
  ggplot(aes(mean.nitrogen,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean 2019 WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 Nitrogen")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# lm 
model_WUE2019_Nitrogen2019 = lm(meanGPP.T ~ Species + mean.nitrogen, data = WUE.2019.Nutrient.2019.data)
model_WUE2019_Nitrogen2019_lmer = lmer(meanGPP.T ~ mean.nitrogen + (1|Species), data = WUE.2019.Nutrient.2019.data)
summary(model_WUE2019_Nitrogen2019)
aov(model_WUE2019_Nitrogen2019)


# 5 
# Mean WUE vs 2019 Leaf N
GPP.T_Nitrogen_2019 = WUE.Nutrient.2019.data %>%
  ggplot(aes(mean.nitrogen,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 Nitrogen")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# lm 
model_WUE_Nitrogen2019 = lm(meanGPP.T ~ Species + mean.nitrogen, data = WUE.Nutrient.2019.data)
model_WUE_Nitrogen2019_lmer = lmer(meanGPP.T ~ mean.nitrogen + (1|Species), data = WUE.Nutrient.2019.data)
summary(model_WUE_Nitrogen2019)
aov(model_WUE_Nitrogen2019)

# 6
# Mean WUE vs 2019 Leaf P
GPP.T_Phosphorus_2019 = WUE.Nutrient.2019.data %>%
  ggplot(aes(mean.phosphorus,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 P (mg ",g^-1," D.M.)")))+ theme_bw() +
  xlim(0.6,1.72) +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# lm 
model_WUE_Phosphorus2019 = lm(meanGPP.T ~ Species + mean.phosphorus, data = WUE.Nutrient.2019.data)
model_WUE_Phosphorus2019_lmer = lmer(meanGPP.T ~ mean.phosphorus + (1|Species), data = WUE.Nutrient.2019.data)
summary(model_WUE_Nitrogen2019)
aov(model_WUE_Nitrogen2019)

# 7
# Mean WUE vs 2019 N/P
GPP.T_NP_2019 = WUE.Nutrient.2019.data %>%
  ggplot(aes(mean.NP.ratio,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 N:P (w/w)")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none") #+ 
  #stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
         #  p.accuracy = 0.001, r.accuracy = 0.001)

# lm 
model_WUE_NP2019 = lm(meanGPP.T ~ Species + mean.NP.ratio, data = WUE.Nutrient.2019.data)
model_WUE_NP2019_lmer = lmer(meanGPP.T ~ mean.NP.ratio + (1|Species), data = WUE.Nutrient.2019.data)
summary(model_WUE_NP2019)
aov(model_WUE_NP2019)


# 8
# Mean WUE vs 2019 K
GPP.T_K_2019 = WUE.Nutrient.2019.data %>%
  ggplot(aes(mean.potassium,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean 2019 WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 K (mg ",g^-1," D.M.)")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none") #+ 
  #stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           #p.accuracy = 0.001, r.accuracy = 0.001)

# lm 
model_WUE_K2019 = lm(meanGPP.T ~ Species + mean.potassium, data = WUE.Nutrient.2019.data)
model_WUE_K2019_lmer = lmer(meanGPP.T ~ mean.potassium + (1|Species), data = WUE.Nutrient.2019.data)
summary(model_WUE_K2019)
aov(model_WUE_K2019)

# 9
# Mean WUE vs 2019 Phosphorus
GPP.T_P_2019 = WUE.Nutrient.2019.data %>%
  ggplot(aes(mean.phosphorus,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean 2019 WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 Phosphorus")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none") + 
  stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
           p.accuracy = 0.001, r.accuracy = 0.001)

# lm 
model_WUE_P2019 = lm(meanGPP.T ~ Species + mean.phosphorus, data = WUE.Nutrient.2019.data)
model_WUE_P2019_lmer = lmer(meanGPP.T ~ mean.phosphorus + (1|Species), data = WUE.Nutrient.2019.data)
summary(model_WUE_P2019)
aov(model_WUE_P2019)

# 10
# Mean WUE vs 2019 Leaf Manganese
GPP.T_Mn_2019 = WUE.Nutrient.2019.data %>%
  ggplot(aes(mean.manganese,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("Mean 2019 WUE (g C ",kg^-1," ",H[2],O,")")), 
       x=expression(paste("Mean 2019 Mn (mg ",g^-1," D.M.)")))+ theme_bw() +
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none") #+ 
  #stat_cor(method = "pearson", label.x.npc = 0.55, label.y.npc = 0.95,
          #p.accuracy = 0.001, r.accuracy = 0.001)
# lm 
model_WUE_Mn2019 = lm(meanGPP.T ~ Species + mean.manganese, data = WUE.Nutrient.2019.data)
model_WUE_Mn2019_lmer = lmer(meanGPP.T ~ mean.manganese + (1|Species), data = WUE.Nutrient.2019.data)
summary(model_WUE_Mn2019)
aov(model_WUE_Mn2019)


# 4 panel image WUE vs nutrients 
Plot_WUE_Nutrients_4_Panel = ggarrange(GPP.T_NP_2019, GPP.T_Phosphorus_2019,
                                       GPP.T_K_2019, GPP.T_Mn_2019,
                                       labels = c("A)", "B)", "C)", "D)"),
                                       font.label = list(size = 12),
                               ncol = 2, nrow = 2)


# carbon validation data - 2019 site mean WUE as a function of c[i]/c[a]
install.packages("ggtext")
library("ggtext")

GPP.T_C13_cica = carbon13validationdata_final %>%
  ggplot(aes(ci.ca,meanGPP.T)) +
  geom_point(alpha=0.8, size=2, aes(color=Species, shape = Species)) +
  labs(y=expression(paste("2019 WUE (g C ",kg^-1," ",H[2],O,")")), x = expression(paste("C"["i"],"/C"["a"])))+ 
  theme_bw() + 
  theme(plot.title = element_text(size = 10, family = "Tahoma"),
        text = element_text(size = 10, family = "Tahoma"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.position="none")

# lm 
model_WUE_C13_cica = lm(meanGPP.T ~ Species + ci.ca, data = carbon13validationdata_final)
model_WUE_C13_cica_lmer = lmer(meanGPP.T ~ ci.ca + (1|Species), data = carbon13validationdata_final)
summary(model_WUE_C13_cica)
aov(model_WUE_C13_cica)




# plot swiz wue plots
plot(mask(mean_Swiz_July_18, forestmask_RP, inverse = TRUE), 
     main = "Mean Summer WUE July 18", breaks = c(1.0, 2.0, 3.0, 4.0),
     axes =  FALSE)

# WUE SWIZ WIDE PLOTS
all_summers_mean_plot = plot(mask(Mean_All_Stacks, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
all_summers_sd_plot = plot(mask(SD_All_Stacks, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
all_summers_median_plot = plot(mask(Median_All_Stacks, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
summer_18_plot = plot(mask(Mean_All_Stacks_18, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
summer_19_plot = plot(mask(Mean_All_Stacks_19, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
summer_20_plot = plot(mask(Mean_All_Stacks_20, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
all_june_plot = plot(mask(Mean_All_Stacks_June, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
all_july_plot = plot(mask(Mean_All_Stacks_July, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))
all_Aug_plot = plot(mask(Mean_All_Stacks_Aug, forestmask_RP, inverse = TRUE), col= brewer.pal(9,"RdYlBu"))

