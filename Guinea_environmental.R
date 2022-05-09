##############################################################################################################################################################################################
################################################ Guinea - extract environmental data from buffered villages #################################################################################
##############################################################################################################################################################################################

library(sp)
library(raster)
library(plyr)
library(SDMTools)
library(rgeos)
library(rgdal)

## Load coastal data
#lfv <- read.csv("C:/Users/iimmkfor/Dropbox (LSoHaTM)/Kim_LSHTM/Back_up_R_files/Guinea/LFV_XYcoordinatesFinal.csv")
lfv2 <- read.csv("C:/Users/iimmkfor/Dropbox (LSoHaTM)/Kim_LSHTM/Back_up_R_files/Guinea/Macenta_village_updated_Jan21.csv")
lfv <- read.csv("C:/Users/iimmkfor/Dropbox (LSoHaTM)/Kim_LSHTM/Back_up_R_files/Guinea/LFV/Lassa_merged.csv")

## Reformat datasets
#lfv <- lfv[c("X_longi", "Y_lat", "Village_name")]
#names(lfv) <- c("x", "y", "vill")
lfv2 <- lfv2[c("Longitude", "Latitude", "Location")]
names(lfv2) <- c("x", "y", "vill")

## X and Y coordinates switched for merged data
colnames(lfv) <- c("result", "sex", "age", "subp", "vill", "code", "study", "site", "y", "x")

## Set coordinates
coordinates(lfv) <- c("x", "y")
proj4string(lfv) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
lfv <- spTransform(lfv, CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"))

## Remove duplicate points and assign code
lfv <- remove.duplicates(lfv)
lfv$vill.name <- lfv$vill
lfv$vill <- paste0("v", seq(1, nrow(lfv), by=1))
lfv <- lfv[c("vill", "vill.name")]

## Set coordinates
coordinates(lfv2) <- c("x", "y")
proj4string(lfv2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
lfv2 <- spTransform(lfv2, CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"))

## Load environmental data 
lc <- raster("E:/GIS/Guinea/Copernicus/ProbaV_2017_DiscreteClass_Macenta_utm.tif")
crs(lc) <- "+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"
tc <- raster("E:/GIS/Guinea/Copernicus/ProbaV_2017_TreeCoverFraction_Macenta_utm.tif")

## Group land classes
lc[lc==114] <- 112
lc[lc==116] <- 112
lc[lc==124] <- 122
lc[lc==126] <- 122

# input vector of distances (in meters)
buffer_distances <- c(500, 1000, 2000, 5000, 10000, 20000)

## Extract land cover around houses
for(i in 1:length(buffer_distances)){
  
  ## Extract data within buffer
  lu <- extract(lc, lfv2, buffer=buffer_distances[i])
  
  ## Summarise by proportions
  sm <- lapply(lu, function(x){prop.table(table(x))})
  
  # convert to data frame
  ev <- data.frame(id = rep(lfv2$vill, lapply(sm, length)),
                   cover = names(unlist(sm)),
                   percent = unlist(sm)
  )
  
  ev$dist <- buffer_distances[i]
  
  # write output
  write.csv(ev, paste0("E:/GIS/Guinea/Macenta_LU_extracted_", buffer_distances[i], ".csv"))
  
}

## Set classification codes
## Land categories:
# 112 = evergreen closed forest, 114 = deciduous closed forest, 116 = closed forest unknown, 122 = evergreen open forest, 124 = deciduous open forest, 126 = open forest unknown, 
# 20 = shrubs, 30 = herbaceous vegetation, 50 = urban/ built up  
cover <- c(112, 122, 20, 30, 50)
cat <- c("cf", "of", "shr", "vg", "bt")
classes <- data.frame(cbind(cover, cat))

## Dataset of house IDS
vl <- data.frame(lfv2$vill)
names(vl) <- "vill"

## Load and recategorise all files
files <- dir("E:/GIS/Guinea", pattern="Macenta_LU")
for(i in 1:length(files)){
  
  ## Read file
  ev <- read.csv(paste0("E:/GIS/Guinea/", files[i]))
  ev$X <- NULL
  
  ## Rename land cover
  ev <- merge(ev, classes, by ="cover")
  dst <- mean(ev$dist)
  ev <- ev[c("id", "percent", "cat")]
  
  ## Merge proportions for all households
  for(j in 1:length(cat)){
    
    # Subset land category
    category <- cat[j]
    temp <- subset(ev, cat==category)
    
    # Rename columns
    names(temp) <- c("vill", paste0(cat[j], dst), "cat")
    temp$cat <- NULL
    temp$vill <- factor(temp$vill)
    
    # Merge to data
    if(exists("ev.pr")){
      ev.pr <- merge(ev.pr, temp, by="vill", all=TRUE)
    }
    if(!exists("ev.pr")){
      ev.pr <- merge(vl, temp, by="vill", all.x=TRUE)
    }
    
    # Replace NAs with 0
    ev.pr[is.na(ev.pr)] <- 0
  }
  
}

write.csv(ev.pr, "E:/GIS/Guinea/Macenta_land_props.csv")

########################## Fragmentation statistics ###################################

df <- lfv2

## Loop through all buffer sizes
for(i in 1:length(buffer_distances)){
  
  ## Loop through each point
  for(j in 1:length(df)){
    
    ## Clip points by buffer radius
    pbuf <- gBuffer(df[j,], width=buffer_distances[i])
    buf <- mask(lc, pbuf)
    
    ## Extract fragstats
    fr <- PatchStat(buf, cellsize = 100, latlon = FALSE)
    fr$vill <- factor(df[j,]$vill)
    fr$dist <- buffer_distances[i]
    
    # Merge to data
    if(exists("fr.df")){
      fr.df <- rbind(fr.df, fr)
    }
    if(!exists("fr.df")){
      fr.df <- fr
    }
    
  }
  
  ## Write csv
  write.csv(fr.df, file = paste0("E:/GIS/Guinea/Macenta_frag_", buffer_distances[i], "m.csv"))
  rm(fr.df)
}

############################## Merge fragmentation stats ################################

## Load and recategorise all files
files <- dir("E:/GIS/Guinea", pattern="Macenta_frag")
for(i in 1:length(files)){
  
  ## Read file
  fr <- read.csv(paste0("E:/GIS/Guinea/", files[i]))
  fr$X <- NULL
  
  ## Rename land cover
  fr <- merge(fr, classes, by.x="patchID", by.y ="cover")
  dst <- mean(fr$dist)
  fr <- fr[c("vill", "perim.area.ratio", "shape.index", "frac.dim.index","cat")]
  
  ## Merge proportions for all households
  for(j in 1:length(cat)){
    
    # Subset land category
    category <- cat[j]
    temp <- subset(fr, cat==category)
    
    # Rename columns
    names(temp) <- c("vill", paste0(cat[j], dst, "pa"), paste0(cat[j], dst, "sh"), paste0(cat[j], dst, "fd"), "cat")
    temp$cat <- NULL
    temp$vill <- factor(temp$vill)
    
    # Merge to data
    if(exists("fr.all")){
      fr.all <- merge(fr.all, temp, by="vill", all=TRUE)
    }
    if(!exists("fr.all")){
      fr.all <- merge(vl, temp, by="vill", all.x=TRUE)
    }
    
    # Replace NAs with 0
    fr.all[is.na(fr.all)] <- 0
  }
  
}

write.csv(fr.all, "E:/GIS/Guinea/Macenta_all_frag_stats.csv")


########################## Mean canopy cover ###################################

## Loop through all buffer sizes
for(i in 1:length(buffer_distances)){
    
    ## Clip points by buffer radius
    tr <- extract(tc, df, buffer=buffer_distances[i], fun=mean, df=TRUE)
    tr <- cbind(tr, df$vill)
    
    ## Rename variables
    names(tr) <- c("ID", paste0("tr", buffer_distances[i]), "vill")
    tr$ID <- NULL
    
    # Merge to data
    if(exists("tr.df")){
      tr.df <- merge(tr.df, tr, by="vill")
    }
    if(!exists("tr.df")){
      tr.df <- tr
    }
    
}

write.csv(tr.df, "E:/GIS/Guinea/Macenta_all_tree_cover_mean.csv")

######################## Merge all datasets together #############################

df <- as.data.frame(df)
df <- merge(df, tr.df, by ="vill")
df <- merge(df, fr.all, by ="vill")
df <- merge(df, ev.pr, by ="vill")

write.csv(df, "E:/GIS/Guinea/Macenta_environmental.csv")

## Remove other datasets
rm(tr.df)
rm(fr.all)
rm(ev.pr)

#######################################################################################################################################################
################################################## Guinea Coastal Areas - 2015 Data ###################################################################
#######################################################################################################################################################

## Load environmental data 
lc <- raster("E:/GIS/Guinea/Copernicus/ProbaV_2015_DiscreteClass_Guinea_UTM.tif")
tc <- raster("E:/GIS/Guinea/Copernicus/ProbaV_2015_TreeCover_Guinea.tif")

## Reproject to UTM
tc <- projectRaster(tc, crs="+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs")

## Group land classes
lc[lc==114] <- 112
lc[lc==116] <- 112
lc[lc==124] <- 122
lc[lc==126] <- 122

# input vector of distances (in meters)
buffer_distances <- c(500, 1000, 2000, 5000, 10000, 20000)

## Extract land cover around houses
for(i in 1:length(buffer_distances)){
  
  ## Extract data within buffer
  lu <- extract(lc, lfv, buffer=buffer_distances[i])
  
  ## Summarise by proportions
  sm <- lapply(lu, function(x){prop.table(table(x))})
  
  # convert to data frame
  ev <- data.frame(id = rep(lfv$vill, lapply(sm, length)),
                   cover = names(unlist(sm)),
                   percent = unlist(sm)
  )
  
  ev$dist <- buffer_distances[i]
  
  # write output
  write.csv(ev, paste0("E:/GIS/Guinea/Coastal_LU_extracted_", buffer_distances[i], ".csv"))
  
}

## Set classification codes
## Land categories:
# 112 = evergreen closed forest, 114 = deciduous closed forest, 116 = closed forest unknown, 122 = evergreen open forest, 124 = deciduous open forest, 126 = open forest unknown, 
# 20 = shrubs, 30 = herbaceous vegetation, 50 = urban/ built up  
cover <- c(112, 122, 20, 30, 50)
cat <- c("cf", "of", "shr", "vg", "bt")
classes <- data.frame(cbind(cover, cat))

## Dataset of house IDS
vl <- data.frame(unique(lfv$vill))
names(vl) <- "vill"
vl$vill <- factor(vl$vill)

## Load and recategorise all files
files <- dir("E:/GIS/Guinea", pattern="Coastal_LU")
for(i in 1:length(files)){
  
  ## Read file
  ev <- read.csv(paste0("E:/GIS/Guinea/", files[i]))
  ev$X <- NULL
  
  ## Rename land cover
  ev <- merge(ev, classes, by ="cover")
  dst <- mean(ev$dist)
  ev <- ev[c("id", "percent", "cat")]
  
  ## Merge proportions for all households
  for(j in 1:length(cat)){
    
    # Subset land category
    category <- cat[j]
    temp <- subset(ev, cat==category)
    
    # Rename columns
    names(temp) <- c("vill", paste0(cat[j], dst), "cat")
    temp$cat <- NULL
    temp$vill <- factor(temp$vill)
    
    # Merge to data
    if(exists("ev.pr")){
      ev.pr <- merge(ev.pr, temp, by="vill", all=TRUE)
    }
    if(!exists("ev.pr")){
      ev.pr <- merge(vl, temp, by="vill", all.x=TRUE)
    }
    
    # Replace NAs with 0
    ev.pr[is.na(ev.pr)] <- 0
  }
  
}

write.csv(ev.pr, "E:/GIS/Guinea/Coastal_land_props.csv")

########################## Fragmentation statistics ###################################

df <- lfv

## Loop through all buffer sizes
for(i in 1:length(buffer_distances)){
  
  ## Loop through each point
  for(j in 1:length(df)){
    
    ## Clip points by buffer radius
    pbuf <- gBuffer(df[j,], width=buffer_distances[i])
    buf <- mask(lc, pbuf)
    
    ## Extract fragstats
    fr <- PatchStat(buf, cellsize = 100, latlon = FALSE)
    fr$vill <- factor(df[j,]$vill)
    fr$dist <- buffer_distances[i]
    
    # Merge to data
    if(exists("fr.df")){
      fr.df <- rbind(fr.df, fr)
    }
    if(!exists("fr.df")){
      fr.df <- fr
    }
    
  }
  
  ## Write csv
  write.csv(fr.df, file = paste0("E:/GIS/Guinea/Coastal_frag_", buffer_distances[i], "m.csv"))
  rm(fr.df)
}

############################## Merge fragmentation stats ################################

## Load and recategorise all files
files <- dir("E:/GIS/Guinea", pattern="Coastal_frag")
for(i in 1:length(files)){
  
  ## Read file
  fr <- read.csv(paste0("E:/GIS/Guinea/", files[i]))
  fr$X <- NULL
  
  ## Rename land cover
  fr <- merge(fr, classes, by.x="patchID", by.y ="cover")
  dst <- mean(fr$dist)
  fr <- fr[c("vill", "perim.area.ratio", "shape.index", "frac.dim.index","cat")]
  
  ## Merge proportions for all households
  for(j in 1:length(cat)){
    
    # Subset land category
    category <- cat[j]
    temp <- subset(fr, cat==category)
    
    # Rename columns
    names(temp) <- c("vill", paste0(cat[j], dst, "pa"), paste0(cat[j], dst, "sh"), paste0(cat[j], dst, "fd"), "cat")
    temp$cat <- NULL
    temp$vill <- factor(temp$vill)
    
    # Merge to data
    if(exists("fr.all")){
      fr.all <- merge(fr.all, temp, by="vill", all.x=TRUE)
    }
    if(!exists("fr.all")){
      fr.all <- merge(vl, temp, by="vill", all.x=TRUE)
    }
    
    # Replace NAs with 0
    fr.all[is.na(fr.all)] <- 0
  }
  
}

write.csv(fr.all, "E:/GIS/Guinea/Coastal_all_frag_stats.csv")


########################## Mean canopy cover ###################################

## Loop through all buffer sizes
for(i in 1:length(buffer_distances)){
  
  ## Clip points by buffer radius
  tr <- extract(tc, df, buffer=buffer_distances[i], fun=mean, df=TRUE)
  tr <- cbind(tr, df$vill)
  
  ## Rename variables
  names(tr) <- c("ID", paste0("tr", buffer_distances[i]), "vill")
  tr$ID <- NULL
  
  # Merge to data
  if(exists("tr.df")){
    tr.df <- merge(tr.df, tr, by="vill")
  }
  if(!exists("tr.df")){
    tr.df <- tr
  }
  
}

write.csv(tr.df, "E:/GIS/Guinea/Coastal_all_tree_cover_mean.csv")

######################## Merge all datasets together #############################

df <- as.data.frame(df)
df <- merge(df, tr.df, by ="vill")
df <- merge(df, fr.all, by ="vill")
df <- merge(df, ev.pr, by ="vill")

write.csv(df, "E:/GIS/Guinea/All_LFV_environmental.csv")
