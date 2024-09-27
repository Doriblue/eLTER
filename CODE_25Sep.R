libs <- c(
  "ncdf4",           # package for netCDF manipulation
  "raster",          # package for raster manipulation
  "ggplot2",         # package for plotting
  "dplyr",           # package for data manipulation
  "basicPlotteR",    # progress bar
  "readxl",          # package for reading Excel files
  "ncdf4.helpers",   # helpers for netCDF manipulation
  "SPEI",            # package for calculating SPEI
  "readr",           # package for reading data
  "tidyr",           # package for data tidying
  "purrr",           # package for functional programming
  "geosphere",       # package for geospatial calculations
  "ncf",             # package for spatial statistics
  "reshape2",        # package for reshaping data
  "sf",
  "cruts"
)

# Install and load libraries
for (pkg in libs) {
  if (!pkg %in% installed.packages()[,"Package"]) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Clean up
rm(pkg)

########FXING NAMES OF SITES##############
renamed<- read_csv("name_fixed.csv")
data<- read_csv("df_masting_PaFs_stValues(fixing site name).csv")
data <- data %>%
  group_by(Latitude, Longitude) %>%
  mutate(Site = first(na.omit(Site)))
colnames(renamed)<- c("Latitude", "Longitude", "Site", "Species","Blank", "SITE")
colnames(data)
# Merge the two dataframes based on the common columns
merged_df <- merge(data, renamed, by = c("Latitude", "Longitude", "Site", "Species"), all.x = TRUE)
data$SITE <- merged_df$SITE
data$Site<- data$SITE
write.csv(data, "df_masting_PaFs_stValues(name_fixed).csv")
Mastree<- read.csv("df_masting_PaFs_stValues(name_fixed).csv")%>%
  filter(Length >= 7) %>%
  group_by(Latitude, Longitude) %>%
  mutate(SITE = ifelse(is.na(SITE), first(na.omit(SITE)), SITE))
# coordinate systems
laea <- "+proj=longlat +ellps=GRS80 +no_defs"
utm <- "+proj=utm +zone=32 +datum=WGS84"
latlong <- "+proj=longlat +ellps=WGS84"
lambert <- "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

###load NUTS###
nuts <-read_sf(dsn = "C:/Users/fpt/Downloads/Vacchiano/eLTER/eLTER/NUTS_2013_60M_SH/NUTS_2013_60M_SH/data", layer =  "NUTS_RG_60M_2013")
nuts <- nuts [nuts$STAT_LEVL_ ==1, "NUTS_ID"]
eu <- read_sf(dsn = "C:/Users/fpt/Downloads/Vacchiano/eLTER/eLTER/ne_10m_admin_0_countries", layer = "ne_10m_admin_0_countries")
eu <- eu[eu$CONTINENT == "Europe", ]
# reproject to UTM 
nuts <- st_transform (nuts, crs = utm)
eu <- st_transform (eu, crs = utm)
# add non-EU countries to NUTS-1
add <- eu [eu$NAME %in% c ("Bosnia and Herz.","Serbia","Ukraine"), "NAME"]
add$NUTS_ID <- c ("BA0", "SR0", "UA0")
add$NAME <- NULL
nuts <- rbind (nuts, add)
# sort alphabetically
nuts <- nuts [order (nuts$NUTS_ID), ]
# calculate centroids 
coord <- st_coordinates(nuts)
nutsGeo <- st_transform(nuts, crs = latlong)
coords <- st_coordinates(st_centroid(nutsGeo))
coords<- as.data.frame(coords)
coords$NUTS_ID<- nutsGeo$NUTS_ID
nutsGeo_spdf <- as(nutsGeo, "Spatial")  # This works because the geometry is a multipolygon
range = c ("1900-01-01","2024-12-31")

##temperature
nc_temp <- nc_open('climatic_data/Tmean_month.nc')
lon <- ncvar_get(nc_temp, "longitude")
lat <- ncvar_get(nc_temp, "latitude", verbose = F)
t <- ncvar_get(nc_temp, "time")
T.array <- ncvar_get(nc_temp, "tg") # store the data in a 3-dimensional array
dim(T.array) 
fill_temp <- ncatt_get(nc_temp, "tg", "_FillValue")
T.array[T.array == fill_temp$value] <- NA
T.slice <- T.array[, , 1] 
dim(T.slice)
r <- raster(t(T.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
plot(r)
writeRaster(r, "TEMPERATURES_NUTS.tif", "GTiff", overwrite=TRUE)
r_brick <- brick(T.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
temp_brick <- flip(t(r_brick), direction='y')
coords$X <- as.numeric(coords$X)
coords$Y <- as.numeric(coords$Y)
dataframe_T <- data.frame(matrix(ncol = length(coords$NUTS_ID), nrow = length(t)))
colnames(dataframe_T) <- coords$NUTS_ID
rownames(dataframe_T) <- TIME
#j<- "DEE"
for (j in coords$NUTS_ID){
  
  #extract data from each point
  pos <- grep(paste0("\\b",j,"\\b"),coords$NUTS_ID)
  lon_pro <- coords$X[pos]
  lat_pro <- coords$Y[pos]
  series_pro <- raster::extract(temp_brick, SpatialPoints(cbind(lon_pro,lat_pro)), method='simple')
  series_pro <- as.vector(series_pro)
  
  dataframe_T[,pos] <- series_pro
  
  # progress(pos, length(coords$NUTS_ID))
  
}


TIME<-nc.get.time.series(f = nc_data,time.dim.name = "time")

dataframe_T <- cbind(TIME = substr(TIME,1,7), dataframe_T) 

write.csv(dataframe_T, "T_NUTS.csv") ##Checked, well referenced

##precipitation
nc_preci <- nc_open('climatic_data/Pcum_month.nc')
lon <- ncvar_get(nc_preci, "longitude")
lat <- ncvar_get(nc_preci, "latitude", verbose = F)
t <- ncvar_get(nc_preci, "time")
P.array <- ncvar_get(nc_preci, "rr") # store the data in a 3-dimensional array
dim(P.array) 
fill_preci <- ncatt_get(nc_preci, "rr", "_FillValue")
P.array[P.array == fill_preci$value] <- NA
P.slice <- P.array[, , 1] 
dim(P.slice)
r <- raster(t(P.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
plot(r)
writeRaster(r, "PRECIPITATION_NUTS.tif", "GTiff", overwrite=TRUE)
r_brick <- brick(P.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
preci_brick <- flip(t(r_brick), direction='y')
coords$X <- as.numeric(coords$X)
coords$Y <- as.numeric(coords$Y)
dataframe_P <- data.frame(matrix(ncol = length(coords$NUTS_ID), nrow = length(t)))
colnames(dataframe_P) <- coords$NUTS_ID
#j<- "DEE"
for (j in coords$NUTS_ID){
  
  #extract data from each point
  pos <- grep(paste0("\\b",j,"\\b"),coords$NUTS_ID)
  lon_pro <- coords$X[pos]
  lat_pro <- coords$Y[pos]
  series_pro <- raster::extract(preci_brick, SpatialPoints(cbind(lon_pro,lat_pro)), method='simple')
  series_pro <- as.vector(series_pro)
  
  dataframe_P[,pos] <- series_pro
  
  # progress(pos, length(coords$NUTS_ID))
  
}


TIME<-nc.get.time.series(f = nc_data,time.dim.name = "time")

dataframe_P <- cbind(TIME = substr(TIME,1,7), dataframe_P) 

write.csv(dataframe_P, "P_NUTS.csv") ##Checked, well referenced

temp <- cruts2poly(ncfile = 'climatic_data/Tmean_reduced.nc', poly = nutsGeo_spdf, timeRange = range, na.rm = TRUE)
preci<- cruts2poly(ncfile = 'climatic_data/Pcum_reduced.nc', poly = nutsGeo_spdf, timeRange = range, na.rm = TRUE)

temp <- as.data.frame(t(temp@data))
preci <- as.data.frame (t (preci@data))

colnames (temp) <- nuts$NUTS_ID
colnames (preci) <- nuts$NUTS_ID

coords <- st_coordinates(st_centroid(nutsGeo))
latitudes <- coords[, 2]  # Latitude is the second column
names(latitudes) <- nutsGeo_spdf$NUTS_ID
###### calculate SPEI#####
pet<- thornthwaite(temp, latitudes, na.rm = T)

balance <- preci-pet
TIME<-nc.get.time.series(f = nc_preci,time.dim.name = "time")
rownames(balance) <- substr(TIME,1,7)

# correct errors in default spei() function
spei_new1 <- eval (parse (text = sub ("-0.35, 0", 
                                      "A=-0.35,B= 0",
                                      deparse (SPEI::spei)), 
                          keep.source = F))

spei_new2 <- eval (parse (text = sub ("s in 1:m", "s in m:1",
                                      deparse (spei_new1)), 
                          keep.source = F))

reassignInPackage("spei", pkgName="SPEI", spei_new2)

# calculate spi and spei
spei3 <- SPEI::spei (balance, scale =3, na.rm =T)
















###########OPEN TEMP###############
  nc_data <- nc_open('climatic_data/Tmean_month.nc')
  
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- ncvar_get(nc_data, "time")
  head(lon) # look at the first few entries in the longitude vector
  T.array <- ncvar_get(nc_data, "tg") # store the data in a 3-dimensional array
  dim(T.array) 
  
  fillvalue <- ncatt_get(nc_data, "tg", "_FillValue")
  fillvalue
  
  nc_close(nc_data)

#CHECK IN THE NETCDF IS GOOD

  T.array[T.array == fillvalue$value] <- NA
  
  T.slice <- T.array[, , 1] 
  
  dim(T.slice)
  
  r <- raster(t(T.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  r <- flip(r, direction='y')
  
  plot(r)
  
  writeRaster(r, "TEMPERATURES.tif", "GTiff", overwrite=TRUE)

  r_brick <- brick(T.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r_brick <- flip(t(r_brick), direction='y')

# EXTRACT ONE DATA TO CHECK
  # example_lon <- 16.67865978
  # example_lat <- 47.75811987
  # example_series <- extract(r_brick, SpatialPoints(cbind(example_lon,example_lat)), method='simple')
  # 
  # example_df <- data.frame(year= seq(from=1, to=length(t), by=1), TEMP=t(example_series))
  # ggplot(data=example_df, aes(x=year, y=TEMP, group=1)) +
  #   geom_line() + # make this a line plot
  #   ggtitle("TEMPERTURE") +     # Set title
  #   theme_bw() # use the black and white theme
  



###### Extract#####
#PiceaAbies
  #PiceaAbies <- read_excel("Pa_coord_fixed.xlsx", 1)
  PiceaAbies <- Mastree%>% filter(Species=="Picea abies")
  PiceaAbies<- PiceaAbies[,c("Site","SITE", "Latitude", "Longitude")] %>% distinct()
  PiceaAbies<- PiceaAbies%>% ungroup %>%
    mutate(PiceaAbies, SITE=paste0(SITE, "_", row_number()))
  PiceaAbies2 <- PiceaAbies
  PiceaAbies[] <- lapply(PiceaAbies, function(x) gsub("\\?", "-", x))
  class(PiceaAbies)
  str(PiceaAbies)
  PiceaAbies$Latitude <- as.numeric(PiceaAbies$Latitude)
  PiceaAbies$Longitude <- as.numeric(PiceaAbies$Longitude)
  
  dataframe_T_Pa <- data.frame(matrix(ncol = length(PiceaAbies$SITE), nrow = length(t)))
  colnames(dataframe_T_Pa) <- PiceaAbies$SITE
  
  for (j in PiceaAbies$SITE){
    
    #extract data from each point
    pos <- grep(paste0("\\b",j,"\\b"),PiceaAbies$SITE)
    lon_pro <- PiceaAbies$Longitude[pos]
    lat_pro <- PiceaAbies$Latitude[pos]
    series_pro <- extract(r_brick, SpatialPoints(cbind(lon_pro,lat_pro)), method='simple')
    series_pro <- as.vector(series_pro)
    
    dataframe_T_Pa[,pos] <- series_pro
    
    # progress(pos, length(PiceaAbies$Site))
    
  }
  
  
  TIME<-nc.get.time.series(f = nc_data,time.dim.name = "time")
  
  dataframe_T_Pa <- cbind(TIME = substr(TIME,1,7), dataframe_T_Pa) 
  
  write.csv(dataframe_T_Pa, "T_Pa.csv") ##Checked, well referenced
  
#Fagus sylvatica 
  #Fagussylvatica <- read_excel("Fs_coord_fixed.xlsx", 1)
  
  Fagussylvatica <- Mastree%>% filter(Species=="Fagus sylvatica")
  Fagussylvatica<- Fagussylvatica[,c("Site","SITE", "Latitude", "Longitude")] %>% distinct()
  Fagussylvatica<- Fagussylvatica%>% ungroup %>%
    mutate(Fagussylvatica, SITE=paste0(SITE, "_", row_number()))
  class(Fagussylvatica)
  str(Fagussylvatica)
  Fagussylvatica$Latitude <- as.numeric(Fagussylvatica$Latitude)
  Fagussylvatica$Longitude <- as.numeric(Fagussylvatica$Longitude)
  
  dataframe_T_Fs <- data.frame(matrix(ncol = length(Fagussylvatica$SITE), nrow = length(t)))
  colnames(dataframe_T_Fs) <- Fagussylvatica$SITE
  
  for (j in Fagussylvatica$SITE){
    
    #extract data from each point
    pos <- grep(paste0("\\b",j,"\\b"),Fagussylvatica$SITE)
    lon_pro <- Fagussylvatica$Longitude[pos]
    lat_pro <- Fagussylvatica$Latitude[pos]
    series_pro <- extract(r_brick, SpatialPoints(cbind(lon_pro,lat_pro)), method='simple')
    series_pro <- as.vector(series_pro)
    
    dataframe_T_Fs[,pos] <- series_pro
    
    # progress(pos, length(Fagussylvatica$Site))
    
  }
  TIME<-nc.get.time.series(f = nc_data,time.dim.name = "time")
  
  dataframe_T_Fs <- cbind(TIME = substr(TIME,1,7), dataframe_T_Fs)
  
  write.csv(dataframe_T_Fs, "T_Fs.csv") ##Checked, well referenced
  
###########OPEN PRECI###############
  nc_data <- nc_open('climatic_data/Pcum_month.nc')
  
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- ncvar_get(nc_data, "time")
  head(lon) # look at the first few entries in the longitude vector
  P.array <- ncvar_get(nc_data, "rr") # store the data in a 3-dimensional array
  dim(P.array) 
  
  fillvalue <- ncatt_get(nc_data, "rr", "_FillValue")
  fillvalue
  
  nc_close(nc_data)
  
  #CHECK IN THE NETCDF IS GOOD
  # P.array[P.array == fillvalue$value] <- NA
  # P.slice <-P.array[, , 1] 
  # dim(P.slice)
  # r <- raster(t(P.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # r <- flip(r, direction='y')
  # plot(r)
  # writeRaster(r, "PRECIPITATION", "GTiff", overwrite=TRUE)
  # 
  
  # EXTRACT ONE DATA TO CHECK
  # r_brick <- brick(P.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # r_brick <- flip(t(r_brick), direction='y')
  # 
  # example_lon <- 16.67865978
  # example_lat <- 47.75811987
  # example_series <- extract(r_brick, SpatialPoints(cbind(example_lon,example_lat)), method='simple')
  # 
  # example_df <- data.frame(year= seq(from=1, to=length(t), by=1), PRECI=t(example_series))
  # ggplot(data=example_df, aes(x=year, y=PRECI, group=1)) +
  #   geom_line() + # make this a line plot
  #   ggtitle("PRECIPITATION") +     # Set title
  #   theme_bw() # use the black and white theme
  
  ###### Extract#####
  #PiceaAbies
  dataframe_P_Pa <- data.frame(matrix(ncol = length(PiceaAbies$SITE), nrow = length(t)))
  colnames(dataframe_P_Pa) <- PiceaAbies$SITE
  
  for (j in PiceaAbies$SITE){
    
    #extract data from each point
    pos <- grep(paste0("\\b",j,"\\b"),PiceaAbies$SITE)
    lon_pro <- PiceaAbies$Longitude[pos]
    lat_pro <- PiceaAbies$Latitude[pos]
    series_pro <- extract(r_brick, SpatialPoints(cbind(lon_pro,lat_pro)), method='simple')
    series_pro <- as.vector(series_pro)
    
    dataframe_P_Pa[,pos] <- series_pro
    
    # progress(pos, length(PiceaAbies$SITE))
    
  }
  
  
  TIME<-nc.get.time.series(f = nc_data,time.dim.name = "time")
  
  dataframe_P_Pa <- cbind(TIME = substr(TIME,1,7), dataframe_P_Pa) 
  
  write.csv(dataframe_P_Pa, "P_Pa.csv") ##Checked, well referenced
  
  #Fagus sylvatica 
  dataframe_P_Fs <- data.frame(matrix(ncol = length(Fagussylvatica$SITE), nrow = length(t)))
  colnames(dataframe_P_Fs) <- Fagussylvatica$SITE
  
  for (j in Fagussylvatica$SITE){
    
    #extract data from each point
    pos <- grep(paste0("\\b",j,"\\b"),Fagussylvatica$SITE)
    lon_pro <- Fagussylvatica$Longitude[pos]
    lat_pro <- Fagussylvatica$Latitude[pos]
    series_pro <- extract(r_brick, SpatialPoints(cbind(lon_pro,lat_pro)), method='simple')
    series_pro <- as.vector(series_pro)
    
    dataframe_P_Fs[,pos] <- series_pro
    
    # progress(pos, length(Fagussylvatica$SITE))
    
  }
  TIME<-nc.get.time.series(f = nc_data,time.dim.name = "time")
  
  dataframe_P_Fs <- cbind(TIME = substr(TIME,1,7), dataframe_P_Fs)
  
  write.csv(dataframe_P_Fs, "P_Fs.csv") ##Checked, well referenced
  
###### calculate SPEI#####
  calculate_pet<- calculate_pet <- function(temp, lat) {
    return(thornthwaite(temp, lat = lat, na.rm = T))
  }
  ###Fagus####
  locations_fs <- Fagussylvatica$SITE
  spei3_fs <- list()
  for (location in locations_fs) {
    precip <- dataframe_P_Fs[[location]]
    temp <- dataframe_T_Fs[[location]]
    lat <- Fagussylvatica[Fagussylvatica$SITE == location, "Latitude"]
    pet <- calculate_pet(temp, lat)
    water_balance <- precip - pet
    spei3 <- spei(water_balance, scale = 3, na.rm = T)
    spei3_fs[[location]] <- spei3
  }
  ###Picea####
  locations_pa <- PiceaAbies$SITE
  spei3_pa <- list()
  for (location in locations_pa) {
    precip <- dataframe_P_Pa[[location]]
    temp <- dataframe_T_Pa[[location]]
    lat <- PiceaAbies[PiceaAbies$SITE == location, "Latitude"]
    pet <- calculate_pet(temp, lat)
    water_balance <- precip - pet
    spei3 <- spei(water_balance, scale = 3, na.rm = T)
    spei3_pa[[location]] <- spei3
  }
  

  #############MASTREE############
  Mastree<- read.csv("df_masting_PaFs_stValues(name_fixed).csv")%>%
    filter(Length >= 7) %>%
    group_by(Latitude, Longitude) %>%
    mutate(SITE = ifelse(is.na(SITE), first(na.omit(SITE)), SITE))
   Mastree_Pa<- Mastree[, c("Year", "Species", "Site","SITE", "Latitude", "Longitude","Start", "End","Length", "Value_st")] %>%
     filter(Species=="Picea abies")%>%
     filter(Year >= 1950)
  Pa <- merge(Mastree_Pa, PiceaAbies, by = c("Latitude", "Longitude", "Site"), all.x = TRUE)
  Pa<- Pa%>% select(-c("Site"))
  Pa<- rename(Pa,Site=SITE)
  Pa_filtered<- Pa %>% 
    filter(!is.na(Latitude) & !is.na(Longitude) & !is.na(Value_st))
  Pa_filtered$Value_st <- as.numeric(Pa_filtered$Value_st)
  
  
  Mastree_Fs<-  Mastree[, c("Year", "Species", "Site","SITE", "Latitude", "Longitude","Start", "End","Length", "Value_st")] %>%
    filter(Species=="Fagus sylvatica")%>%
    filter(Year >= 1950)
  Fs <- merge(Mastree_Fs, Fagussylvatica, by = c("Latitude", "Longitude", "Site"), all.x = TRUE)
  Fs<- Fs%>% select(-c("Site"))
  Fs<- rename(Fs,Site=SITE)
  Fs_filtered <- Fs %>%
    filter(!is.na(Latitude) & !is.na(Longitude) & !is.na(Value_st))
  Fs_filtered$Value_st <- as.numeric(Fs_filtered$Value_st)
###MANTEL test-----------
##FS----
  Fs_filtered_agg <- aggregate(Value_st ~ Year + Site, data = Fs_filtered, FUN = mean, na.rm = TRUE)
  Fs_wide <- Fs_filtered_agg %>%
    pivot_wider(names_from = Site, values_from = Value_st)
  print(str(Fs_wide)) 
  cor_data <- Fs_wide[,-1]  
  if (ncol(cor_data) > 0) {
    spear.o <- cor(cor_data, use = "pairwise.complete.obs", method = "spearman")
    name.nuts <- colnames(cor_data) 
    rownames(spear.o) <- name.nuts
    colnames(spear.o) <- name.nuts
      for (i in name.nuts) {
      spear.o[i, i] <- NA
    }}
    
  include <- names(which(apply(spear.o, 2, sum, na.rm = TRUE) != 0))
  spear.o <- spear.o[rownames(spear.o) %in% include, colnames(spear.o) %in% include]
  
  d.o <- as.dist((-spear.o + 1) / 2 * 100) #dissimilarituy distance metrix of masting between sites
  sites_coords <- Fs_filtered[Fs_filtered$Site %in% rownames(spear.o), ]
  sites_coords<- sites_coords[,c(1,2,5)]%>% distinct(Site, .keep_all = TRUE)
  coord <- as.matrix(sites_coords[, c("Longitude", "Latitude")])
  
  # Distance matrix
  dgeo <- dist(coord)
  
  # Step 6: Plot similarity vs distance
  dgeo.m <- as.matrix(dist(coord, diag = TRUE, upper = TRUE))
  controprova <- data.frame("dgeo" = as.vector(dgeo.m*100),
                            "similarity" = as.vector(spear.o))
  controprova <- controprova[order(controprova$dgeo, decreasing = TRUE), ]
  
  # Scatter plot
  plot(controprova$dgeo, controprova$similarity, 
       col = "dark grey", pch = 16, cex = 0.5, 
       ylab = "Similarity", xlab = "Distance [km]", las = 1)
    abline(h = 0, col = "grey")
    loess_fit <- loess(similarity ~ dgeo, data = controprova)
    distance_seq <- seq(min(controprova$dgeo), max(controprova$dgeo), length.out = 100)
    predicted_values <- predict(loess_fit, newdata = data.frame(dgeo = distance_seq))
    lines(distance_seq, predicted_values, col = "blue", lwd = 3)
  
  # Step 7: Mantel test
  mantel_test <- ncf::mantel.test(M1 = as.matrix(dgeo*100), 
                                  M2 = as.matrix(d.o), 
                                  resamp = 2000, quiet = TRUE)[1:2]  # $correlation[1] 0.006084475 $p[1] 0.3883058
  
  # Step 8: Mantel correlogram
  plot_mantel <- mantel.correlog(as.matrix(dgeo*100), 
                                 as.matrix(as.dist(spear.o)), 
                                 increment = 400000, 
                                 resamp = 2000, 
                                 quiet = TRUE)
  
  # Plot correlogram
  plot(plot_mantel$mean.of.class / 1000, plot_mantel$correlation, 
       xlim = c(0, 1), 
       xlab = "Distance [km]", 
       ylab = "Spearman's correlation", type = "b", 
       pch = ifelse((is.na(plot_mantel$p) | plot_mantel$p > 0.05), 1, 16), las = 1)
  abline(h = 0, col = "grey")
  
##PA----
  
  
# redo with distance along latitude ad longitude only
# by fixing the other coordinate to mean value
coordLat <- coord
coordLat [, 1] <- mean (coordLat [, 1])
dgeoLat <- dist (coordLat)

coordLong <- coord
coordLong [, 2] <- mean (coordLong [, 2])
dgeoLong <- dist (coordLong)

ncf::mantel.test (M1 = as.matrix (dgeoLat), M2 = as.matrix (d.o), resamp = 2000, quiet = T) [1:2] #$correlation[1] 0.03277124 $p[1] 0.07896052
ncf::mantel.test (M1 = as.matrix (dgeoLong),M2 = as.matrix (d.o), resamp = 2000, quiet = T) [1:2] # #$correlation[1] -0.003403767 $p[1] 0.4387806

plot.lat <- mantel.correlog (as.matrix (dgeoLat), as.matrix (as.dist (spear.o)), increment = 400000, resamp = 2000, quiet=T)
plot.long <- mantel.correlog (as.matrix (dgeoLong),as.matrix (as.dist (spear.o)), increment = 400000, resamp = 2000, quiet=T)

par (mfcol = c (1, 2))
plot (plot.lat$mean.of.class / 1000, plot.lat$correlation, 
      ylim = c (-0.1, 0.5), xlim = c (0, 2500),
      xlab = "distance over latitude [km]",
      ylab = "Spearman's correlation", type = "b", 
      pch = ifelse ((is.na(plot.lat$p) | plot.lat$p>0.05), 1, 16), las =1)
abline (h = 0, col = "grey")
plot (plot.long$mean.of.class / 1000, plot.long$correlation, 
      ylim = c (-0.1, 0.5), xlim = c (0, 2500),
      xlab = "distance over longitude [km]",
      ylab = "Spearman's correlation", type = "b", 
      pch = ifelse ((is.na(plot.long$p) | plot.long$p>0.05), 1, 16), las = 1)
abline (h =0, col="grey")
par (mfcol = c (1, 1))


### hierarchical clustering ####

# fit rho-distance model 
dgeo.mLong <- as.matrix (dist (coordLong, diag = T, upper = T))
dgeo.mLat <-  as.matrix (dist (coordLat,  diag = T, upper = T))
df2 <- data.frame ("dgeoLong" = as.vector (dgeo.mLong),
                   "dgeoLat" = as.vector (dgeo.mLat),
                   "dissimilarity" = as.vector (as.matrix(d.o)))

model <- lm (dissimilarity ~ dgeoLong + dgeoLat, data = df2)

df2$intercept <- model$coefficients[ 1]
df2$bLong <- model$coefficients [2]
df2$bLat <- model$coefficients [3]
df2$fitted <- df2$intercept + df2$bLong * df2$dgeoLong + df2$bLat * df2$dgeoLat

# calculate ratio of NAs
# d.o.full <- as.vector ((-spear.o + 1) / 2 * 100)
n.na <- length (d.o [is.na (d.o)])
ratio <- n.na / length (d.o) # 4.9%

# apply fitted values to NAs in similarity matrix
df2 <- df2 [which (is.na (df2$dissimilarity)), ]
df2 <- df2 [df2$dgeoLat > 0, ]
df2$check = 1

for (i in 1: (nrow (df2) - 1)){
  for (j in (i+1) : nrow (df2)){
    if (df2$fitted [i] == df2$fitted [j]) {df2$check [j] <- 0}
  }
}

df2 <- df2 [df2$check == 1, ]
d.o.fitted <- d.o
d.o.fitted [is.na (d.o)] <- df2$fitted

# clustering with Ward D2 distance
cluster.o <- hclust (d.o.fitted, method = "ward.D2")

# optimal number of clusters
level.o <- NbClust (diss = d.o.fitted, distance = NULL, index="mcclain",
                    min.nc = 3, max.nc = 15, method = "ward.D2")
level <- level.o$Best.nc [1] #3

# bootstrap clusters for stability
cboot.hclust.o <- clusterboot (d.o.fitted, clustermethod = hclustCBI,
                               distances = T, method = "ward.D2", k = level)
cboot.hclust.o$bootmean # 0.7188821 0.5437268 0.6478351

# plot dendrograms
par (mfcol = c (1, 1))
plot (cluster.o, cex = 0.7, ylab = "Suzuki dissimilarity")
rect.hclust (cluster.o, k = level)

# define clusters
mycl.o <- data.frame ("cluster"=  cutree (cluster.o, k = level))
mycl.o$id <- names (cutree (cluster.o, k = level))

# map clusters
final.plot.o <- merge (nutsFasyMap, mycl.o, by = "id", all.x = TRUE)

ggplot () +
  geom_polygon (data = eu,
                aes(x = long, y = lat, group = group),
                size = 0) + 
  xlim (2500000, 5900000) +
  ylim (1500000, 4500000) +
  coord_equal () +
  
  geom_polygon (data = final.plot.o, 
                aes (x = long, y = lat, group = group, fill = factor (cluster)), 
                color = NA,size = 0) + 
  theme (legend.position = "none") 


# data_spei_fs <- list()
# for (site in names(spei3_fs)) {
#   data_spei_fs[[site]] <- spei3_fs[[site]][["fitted"]]
# }
# spei_fs_df <- do.call(cbind, data_spei_fs)
# TIME_Fs<- dataframe_P_Fs$TIME
# spei_fs_df <- data.frame(TIME_Fs, spei_fs_df)
# 
# colnames(spei_fs_df) <- c("TIME", names(data_spei_fs))
# 
# head(spei_fs_df)
# 
# 
# spei_fs_df$TIME <- as.Date(paste0(spei_fs_df$TIME, "-01"), format = "%Y-%m-%d")
# 
# extract_spei3_series <- function(start_year, end_year, site) {
#   start_date <- as.Date(paste0(start_year, "-01-01"))
#   end_date <- as.Date(paste0(end_year, "-12-31"))
#   spei3_series <- spei_fs_df %>%
#     filter(TIME >= start_date & TIME <= end_date) %>%
#     select(TIME, all_of(site))
#   return(spei3_series[[site]])
# }
# 
# # Extract spei3 series and add to Fs
# Fs <- Fs %>%
#   rowwise() %>%
#   mutate(spei3_series = list(extract_spei3_series(Start, End, Site)))
# 
# # Unnest the spei3 series
# Fs_unnested <- Fs %>%
#   unnest(cols = c(spei3_series)) %>%
# filter(!is.na(spei3_series) & !is.na(Value_st))
# 
# # Create distance matrices
# spei3_dist <- dist(Fs_unnested$spei3_series, method = "euclidean")
# masting_dist <- dist(Fs_unnested$Value_st, method = "euclidean")
# spei3_dist_matrix <- as.matrix(spei3_dist)
# masting_dist_matrix <- as.matrix(masting_dist)
# mantel_result <- vegan::mantel(spei3_dist_matrix, masting_dist_matrix, permutations = 2000)
