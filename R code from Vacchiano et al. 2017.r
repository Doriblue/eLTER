### Spatio-temporal patterns of beech masting####
# Code by Giorgio Vacchiano, Andy Hacket-Pain# 
# Version 2017-02-12 # 
# Contact info: gvacchiano@gmail.com #

### Cleaning and libraries ####
rm (list = ls () )

# Required libraries
lib <- c ("sp", "quantmod", "ncf", "ggplot2", "cluster", "gplots",
          "fpc", "gdata", "NbClust", "tseries", "repmis", "RCurl","rgdal",
          "xlsx", "gnumeric", "multitaper", "psd", "dplR", "gtools", "raster",
          "rgeos", "phytools", "grid", "SpatialPack", "GeneNet", "ordinal",
          "ncdf4", "chron", "cruts", "SPEI", "R.utils", "MASS", "nnet", "rms",
          "vioplot","gam")

# Install and load libraries
for (i in lib) {
  if (i %in% installed.packages () [,"Package"] == F) {
    install.packages (i)
  }
  library (i, character.only = TRUE)
}

rm(i)

### Spatial objects loading ####

# data URLs 
url.nuts <- "http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/NUTS_2013_60M_SH.zip"
url.eu <- "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip"
url.g16a <- "http://ies-ows.jrc.ec.europa.eu/efdac/download/Atlas/Fagus-sylvatica_rpp.zip"
url.pre <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.08/cruts.2406270035.v4.08/pre/cru_ts4.08.1901.2023.pre.dat.nc.gz"
url.tmn <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.08/cruts.2406270035.v4.08/tmn/cru_ts4.08.1901.2023.tmn.dat.nc.gz"
url.tmp <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.08/cruts.2406270035.v4.08/tmp/cru_ts4.08.1901.2023.tmp.dat.nc.gz"
url.tmx <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.08/cruts.2406270035.v4.08/tmx/cru_ts4.08.1901.2023.tmx.dat.nc.gz"
url.pet <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.08/cruts.2406270035.v4.08/pet/cru_ts4.08.1901.2023.pet.dat.nc.gz"

# coordinate systems
laea <- "+proj=longlat +ellps=GRS80 +no_defs"
utm <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
latlong <- "+proj=longlat +ellps=WGS84"
lambert <- "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# load NUTS-1
dir <- paste0 (tempdir (),runif (1, 0, 9))
temp <- tempfile (fileext = ".zip")
download.file (url.nuts, temp)
unzip (temp, exdir = dir)
nuts <- readOGR (dsn = paste0 (dir, "/NUTS_2013_60M_SH/data"), 
                     layer = "NUTS_RG_60M_2013", verbose=FALSE)

nuts <- st_read(dsn = paste0 (dir, "/NUTS_2013_60M_SH/data"), 
                layer = "NUTS_RG_60M_2013")

nuts <- nuts [nuts$STAT_LEVL_ ==1, "NUTS_ID"]

# load countries
dir <- paste0 (tempdir (),runif (1, 0, 9))
temp <- tempfile (fileext = ".zip")
download.file (url.eu, temp)
unzip (temp, exdir = dir)
f <- list.files (path="C:/Users/fpt/Downloads/Vacchiano/eLTER/eLTER/ne_10m_admin_0_countries", ".shp", full.names = T)
eu <- maptools::readShapePoly (f, proj4string = CRS (latlong))
eu <- eu [eu$CONTINENT == "Europe",]

# load  beech distribution map (20 Mb)
dir <- paste0 (tempdir (),runif (1, 0, 9))
temp <- tempfile (fileext = ".zip")
download.file (url.g16a, temp, mode="wb")
unzip (temp, exdir = dir)
g16a <- raster(list.files (dir, ".tif", full.names = T))

# threshold beech map 
g16a [g16a <0.05] <- NA
g16a [g16a>=0.05] <- 1

# convert to Polygon and dissolve
# WARNING VERY LONG (and crash-prone) OPERATIONS
euforgen <- rasterToPolygons (g16a)
euforgen <- gUnaryUnion (euforgen, id=euforgen@data$DN)

# reproject to UTM 
nuts <- spTransform (nuts, CRS = utm)
eu <- spTransform (eu, CRS = utm)
euforgen <- spTransform (euforgen, CRS = utm)

# add non-EU countries to NUTS-1
add <- eu [eu$NAME %in% c ("Bosnia and Herz.","Serbia","Ukraine"), "NAME"]
add$NUTS_ID <- c ("BA0", "SR0", "UA0")
add$NAME <- NULL
nuts <- rbind (nuts, add)

# sort alphabetically
nuts <- nuts [order (nuts$NUTS_ID), ]

# calculate centroids 
coord <- getSpPPolygonsLabptSlots (nuts)

### Climate data extraction ####

gz.tmn <- tempfile (pattern = "tmn", fileext = ".gz")
download.file (url.tmn, gz.tmn); gunzip (gz.tmn) #156Mb
link.tmn <- list.files (dirname (gz.tmn), pattern = "tmn", full.names =T) [1] 

gz.tmx <- tempfile (pattern = "tmx",fileext = ".gz")
download.file (url.tmx, gz.tmx); gunzip (gz.tmx) #156Mb
link.tmx <- list.files (dirname (gz.tmx), pattern = "tmx", full.names =T) [1] 

gz.tmp <- tempfile (pattern = "tmp",fileext = ".gz")
download.file (url.tmp, gz.tmp); gunzip (gz.tmp) #156Mb
link.tmp <- list.files (dirname (gz.tmp), pattern = "tmp", full.names =T) [1] 

gz.pre <- tempfile (pattern = "pre",fileext = ".gz")
download.file (url.pre, gz.pre); gunzip (gz.pre) #199Mb
link.pre <- list.files (dirname (gz.pre), pattern = "pre", full.names =T) [1] 

gz.pet <- tempfile (pattern = "pet",fileext = ".gz")
download.file (url.pet, gz.pet); gunzip (gz.pet) #78Mb
link.pet <- list.files (dirname (gz.pet), pattern = "pet", full.names =T) [1] 

nutsGeo <- spTransform (nuts, CRS = latlong)
range = c ("1900-01-01","2016-12-31")

tmn <- cruts2poly (link.tmn, nutsGeo, timeRange = range, na.rm = T)
tmx <- cruts2poly (link.tmx, nutsGeo, timeRange = range, na.rm = T)
tmp <- cruts2poly (link.tmp, nutsGeo, timeRange = range, na.rm = T)
pre <- cruts2poly (link.pre, nutsGeo, timeRange = range, na.rm = T)
pet <- cruts2poly (link.pet, nutsGeo, timeRange = range, na.rm = T)

tmn <- as.data.frame (t (tmn@data))
tmx <- as.data.frame (t (tmx@data))
tmp <- as.data.frame (t (tmp@data))
pre <- as.data.frame (t (pre@data))
pet <- as.data.frame (t (pet@data))

colnames (tmn) <- nuts$NUTS_ID
colnames (tmx) <- nuts$NUTS_ID
colnames (tmp) <- nuts$NUTS_ID
colnames (pre) <- nuts$NUTS_ID
colnames (pet) <- nuts$NUTS_ID

balance <- pre-pet

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
spi3 <- as.data.frame (SPEI::spi (pre, scale =3, na.rm =T, fit = "pp-pwm")$fitted)
spei3 <- as.data.frame (SPEI::spei (balance, scale =3, na.rm =T)$fitted)

tmn$month   <- rep (1:12, nrow (tmn)/12)
tmx$month   <- rep (1:12, nrow (tmn)/12)
tmp$month   <- rep (1:12, nrow (tmn)/12)
pre$month   <- rep (1:12, nrow (tmn)/12)
pet$month   <- rep (1:12, nrow (tmn)/12)
spi3$month  <- rep (1:12, nrow (tmn)/12)
spei3$month <- rep (1:12, nrow (tmn)/12)

tmn$year   <- substr (rownames (tmn), 2, 5)
tmp$year   <- substr (rownames (tmn), 2, 5)
tmx$year   <- substr (rownames (tmn), 2, 5)
pre$year   <- substr (rownames (tmn), 2, 5)
pet$year   <- substr (rownames (tmn), 2, 5)
spi3$year  <- substr (rownames (tmn), 2, 5)
spei3$year <- substr (rownames (tmn), 2, 5)

rownames (tmn) <- NULL
rownames (tmx) <- NULL
rownames (tmp) <- NULL
rownames (pre) <- NULL
rownames (pet) <- NULL
rownames (spi3) <- NULL
rownames (spei3) <- NULL

# create and bind NA rows for 2015-16 (no climate data)
clim.na.col <- as.data.frame (matrix (rep (NA, (length (nuts) *24)),
                                      ncol = length (nuts)))
clim.na <- cbind (rep (NA, 24), 1:24, clim.na.col)
                    
names (clim.na) <- names (tmx)
clim.2015 <- clim.na
clim.2015$year <- c(rep(2015,12),rep(2016,12))
clim.2015$month <- c(1:12,1:12)

tmn   <- rbind (tmn, clim.2015)
tmx   <- rbind (tmx, clim.2015)
tmp   <- rbind (tmp, clim.2015)
pre   <- rbind (pre, clim.2015)
pet   <- rbind (pet, clim.2015)
spi3  <- rbind (spi3, clim.2015)
spei3 <- rbind (spei3, clim.2015)

### Masting data - loading ####

# data URL
url.data <- "http://onlinelibrary.wiley.com/store/10.1002/ecy.1785/asset/supinfo/ecy1785-sup-0002-DataS1.zip?v=1&s=2491b8cc559d5ec909f96dfc5a91397b1d7e9683"
dir <- paste0 (tempdir (),runif (1, 0, 9))
temp <- tempfile (fileext = ".zip")
download.file (url.data, temp, method="wget", extra="--no-check-certificate")
unzip (temp, exdir = dir)

# read MASTREE data
df <- read.csv (paste0 (dir, "/MASTREE_2016.11.csv"))

# data filtering and transformation
df <- df [df$AccessionDate != "September_2016", ] # filter by accession date (june_2016)
df <- df [df$Species == "FASY", ] # filter by species (beech)

# exclude flower and pollen data
df <- df [df$Proxy =="seed" | df$Proxy =="fruit" | df$Proxy =="cone", ]

# list unique NUTS-1
nutsDef <- unique (df$Nuts1) 
nutsDef <- nutsDef [order (nutsDef)]

# remove territorial units larger than NUTS-1
remove <- nutsDef [substr (nutsDef,3,3)=="#"] 
df <- df[! df$Nuts1 %in% remove, ]

# list unique NUTS-1
nutsDef <- unique (df$Nuts1) 
nutsDef <- nutsDef [order (nutsDef)]
coordDef<-getSpPPolygonsLabptSlots (nuts [nuts$NUTS_ID %in% nutsDef,])

### correlation between SPEI and summer Tmax ####
spei.jj <- aggregate (spei3 [spei3$month %in% c (6,7), 1:(ncol(spei3)-2)],
                   by = list (spei3$year [spei3$month %in% c (6,7)]), 
                   FUN = mean)
spei.jj$year <- unique (spei3$year)

tmax.jj <- aggregate (tmx [tmx$month %in% c (6,7), 1:(ncol(tmx)-2)],
                      by = list (tmx$year [tmx$month %in% c (6,7)]), 
                      FUN = mean)
tmax.jj$year <- unique (tmx$year)

spei.jj <- spei.jj [, colnames (spei.jj) %in% nutsDef]
tmax.jj <- tmax.jj [, colnames (tmax.jj) %in% nutsDef]

cor.jj <- data.frame ("id"=nutsDef, "r" = NA)
for (i in 1:ncol (spei.jj)) {
  cor.jj$r[i] <- cor (spei.jj [, i], tmax.jj [, i], use = "complete.obs")
}

# plot intra-NUTS drought-Tmax correlation
eu <- fortify (eu)
nutsFasyMap <- fortify (nuts, region = "NUTS_ID")
merge.shp.coef <- merge (nutsFasyMap, cor.jj, by = "id")
final.plot <- merge.shp.coef [order (merge.shp.coef$order), ] 

ggplot () +
  geom_polygon (data = eu, # background admin boundaries
                aes (x = long, y = lat, group = group),
                size = 0.25) + 
  xlim (2500000, 6500000) +
  ylim (1500000, 4500000) +
  geom_polygon (data = final.plot, 
                aes (x = long, y = lat, group = group, fill = r), 
                color = "black", size = 0.1) +
  scale_fill_gradient2 (guide = "legend", high = "dark red") 

# scatterplot
spei.jj[spei.jj>5]<-NA

plot (seq (min (spei.jj, na.rm=T), max (spei.jj, na.rm=T), length.out = 10), 
      seq (min (tmax.jj, na.rm=T), max (tmax.jj, na.rm=T), length.out = 10), 
      type="n", xlab="SPEI3 June-July", ylab="Tmax June-July (°C)")

for (i in 1:ncol (spei.jj)) {
  lines (loess.smooth (spei.jj [, i], tmax.jj [, i]), col=i)
}

tmax.ave <- colMeans(tmax.jj, na.rm=T)
spei.jj<- spei.jj[,order(tmax.ave)]
tmax.jj<- tmax.jj[,order(tmax.ave)]

#tiff (filename = "drought2.tiff", width = 1600, height = 1600)
par (mfrow=c(7,9))
for (i in 1:ncol(spei.jj)){
  scatter.smooth (spei.jj [, i], tmax.jj [, i], 
                  pch=16, xlim=c(-2.5,+2.5), ylim =c(15,30),
                  xlab="SPEI JJ", ylab = "Tmax JJ", main = names(spei.jj[i]))
}
dev.off()
par(mfcol=c(1,1))

### climate homogeneity of NUTS1 ####
nutsFasy <- nuts [nuts$NUTS_ID %in% nutsDef,]
polyH=nutsFasy
spTransform(polyH, CRS("+init=epsg:4326"))

#ncfile=link.tmn; outfile="tmnCor.tiff"
#ncfile=link.tmx; outfile="tmxCor.tiff"
#ncfile=link.tmp; outfile="tmpCor.tiff"
#ncfile=link.pre; outfile="preCor.tiff"
#ncfile=link.pet; outfile="petCor.tiff"

br <- cruts2raster(ncfile = ncfile, timeRange = range, 
                   poly = polyH, offset = "1900-01-01" , type = "brick")

corm=data.frame("id"=nutsDef,"y"=NA)
for (i in 1:length(polyH)){
  ext<-raster::extract(br, polyH[i,], df=T)
  ext <- t(as.matrix (ext))
  corm$y[i] <- mean(cor (ext, use="complete.obs") [lower.tri(cor (ext),diag=F)])
}

# plot intra-NUTS climate correlation
merge.shp.coef <- merge (nutsFasyMap, corm, by = "id")
final.plot <- merge.shp.coef [order (merge.shp.coef$order), ] 

# tiff (filename = outfile, width = 800, height = 600)

ggplot () +
  geom_polygon (data = eu, # background admin boundaries
                aes (x = long, y = lat, group = group),
                size = 0.25) + 
  xlim (2500000, 6500000) +
  ylim (1500000, 4500000) +
  geom_polygon (data = final.plot, 
                aes (x = long, y = lat, group = group, fill = y), 
                color = "black", size = 0.1) +
  scale_fill_gradient2 (guide = "legend", high = "dark red") 

#dev.off()

### cluster analysis of climate ####

climL = list (tmp,tmn,tmx,pre,spi3,spei3)
names(climL) <- c("tmp","tmn","tmx","pre","spi3","spei3")

L=1
L=2
L=3
L=4
L=5
L=6

climdf <- climL[[L]]

climdf <- climdf[,names(climdf) %in% nutsDef]
climdf$month <- NULL
climdf$year <- NULL

spear.L <- cor (climdf, use = "pairwise.complete.obs")
for (i in nutsDef) {spear.L [i, i] <- NA}

include <- names (which (apply (spear.L, 2, sum, na.rm = T) != 0))
spear.L <- spear.L [rownames (spear.L) %in% include, colnames (spear.L) %in% include]

d.L.fitted <- as.dist ((-spear.L + 1) / 2 * 100)

cluster.L <- hclust (d.L.fitted, method = "ward.D2")
level.L <-3

mycl.L <- data.frame ("cluster"=  cutree (cluster.L, k = level.L))
mycl.L$id <- names (cutree (cluster.L, k = level.L))
mycl.L[! mycl.L$id %in% nutsDef] <- NA

final.plot.L <- merge (nutsFasyMap, mycl.L, by = "id", all.x = TRUE)

outfile=paste0(names(climL)[L],".cluster.tiff")
tiff (filename = outfile, width = 800, height = 600)

ggplot () +
  geom_polygon (data = eu,
                aes(x = long, y = lat, group = group),
                size = 0) + 
  xlim (2500000, 5900000) +
  ylim (1500000, 4500000) +
  coord_equal () +
  
  geom_polygon (data = final.plot.L, 
                aes (x = long, y = lat, group = group, fill = factor (cluster)), 
                color = "black",size = 0.1) + 
  theme (legend.position = "none") 
dev.off()

### Descriptive statistics of individual series ####

# select columns for Intra-Nuts correlation
head(df)
df.series <- df [, c ("ORDmast", "Yr", "ID", "Length", "Nuts1")]

head(df)
# paste Nuts-1 and ID
df.series$ID2 <- paste (df.series$Nuts1, df.series$ID, sep = "_")

# delete series with length <X years
X = 7
df.series7 <- df.series [df.series$Length >= X, ]

# reshape to wide format
series7 <- reshape (df.series7 [, c("Yr", "ORDmast", "ID2")],
                   direction= "wide",
                   idvar= "Yr", 
                   timevar= "ID2",
                   v.names= "ORDmast")

series <- reshape (df.series [, c("Yr", "ORDmast", "ID2")],
                    direction= "wide",
                    idvar= "Yr", 
                    timevar= "ID2",
                    v.names= "ORDmast")

# sort by year
series7 <- series7 [order (series7$Yr), ]
series <- series [order (series$Yr), ]

# relative frequencies of masting classes
length <- apply (series7 [, -1], 2, function (x) length (na.omit (x)))
n1 <- apply (series7 [, -1], 2, function (x) length (na.omit (x [x==1])))
n2 <- apply (series7 [, -1], 2, function (x) length (na.omit (x [x==2])))
n3 <- apply (series7 [, -1], 2, function (x) length (na.omit (x [x==3])))
n4 <- apply (series7 [, -1], 2, function (x) length (na.omit (x [x==4])))
n5 <- apply (series7 [, -1], 2, function (x) length (na.omit (x [x==5])))

desc <- c (mean (n1/length),
           mean (n2/length),
           mean (n3/length),
           mean (n4/length),
           mean (n5/length)) #0.38 0.25 0.18 0.08 0.10

# ### creating a Skeleton Plot of all individual series ####
# commented out to shorten runtime
# 
# # truncate to required time period
# series.trun <- series [series$Yr >= 1901,]
# 
# # replace NAs with zero (for plotting, they are just plotted white)
# series.trun [is.na (series.trun)] <- 0
# 
# # sort by series ID
# series.trun <- series.trun [,2:ncol (series.trun)][, order (
#   colnames (series.trun [,2:ncol(series.trun)]))]
# 
# # plotting
# png ("SKELETON PLOT.png", width=20, height=56, units="cm", res=800)
# 
# col0<-  "white"
# col20<- "dark orange"
# col1<- "dark green"
# rbPal <- colorRampPalette(c(col0, col20, col1))
# Col <- rbPal(6)
# 
# data<- series.trun
# m <- nrow (data)   #number of NUTs regions
# n <- ncol (data)
# x <- seq (m)
# y <- seq (n)
# 
# shortnames <- sapply (colnames (data), strsplit, "FASY")
# ylabs <- c (); for (i in 1:length (shortnames)) {ylabs [i] <- shortnames[[i]][2]}
# xlabs<- 1901:2016
# 
# par(mar=c(3,4,.4,4)+0.1, mgp=c(1.1,0.2,0), tck=0.02, lwd=1, cex=1, cex.lab=0.5, cex.axis=0.8,  pch=16)
# par(xaxs="i"); par(yaxs="i")
# z <-data.matrix (data) #rownames.force =NA)
# image (x,y,z, col= Col, axes=FALSE, ylab="", xlab="")
# 
# mtext ("YEAR",1, line=1, cex=0.7)
# axis (2, at=c(1:n), tcl=0, labels=ylabs, las=1, cex.axis=0.05)
# axis (4, at=c(1:n), tcl=0, labels="", las=1, cex.axis=0.38)
# axis (1, at=c(1:m), tcl=0, labels=xlabs, las=2, cex.axis=0.3)
# #grid (nx = m, ny=n, col="grey80", lty="solid", lwd=1)
# box (col="grey80", lwd=1)
# 
# dev.off()

### Intra-Nuts correlation ####

# extract unique Nuts1 IDs
series.nuts <- unlist (lapply 
                       (strsplit (
                         sub ("ORDmast.", "", colnames (series) [-1]),
                         "_", fixed = T, useBytes = T),
                       function (l) l [[1]]
                       )
                       )

# prepare vectors for correlation and number of observation pairs
series.nuts.cor<- rep (NA, length (unique (series.nuts)))
series.nuts.n <- rep (NA, length (unique (series.nuts))) 

# calculate intra-nuts correlations
for (i in 1:length (series.nuts.cor)){
  selected <- which (series.nuts == unique (series.nuts) [i])
  if (length (selected) > 1) {
    series.test <- series [,-1] [, selected]
    year <- series$Yr
    cor.matrix <- c ()
    np <- c ()
    for (j in 1:ncol (series.test)){
      for (k in (1:ncol (series.test)) [-c (1:j)]){
        np <- c (np, nrow (na.omit (series.test[, c (j, k)])))
        corJK <- na.omit (cbind (year, series.test [,c (j, k)]))
        # filter pairs with >= X years of measurement
        if (tail (np,1) <X) {kk <- NA} 
        if (tail (np,1) >=X) {
          if (sd (corJK [, 2]) == 0 | sd (corJK [, 3]) == 0) {kk <- NA} 
          else {
            # adjusted degrees fo freedom (Dutilleuil 1963)
            kk <- modified.ttest  (corJK[, 2], corJK[, 3],
                                   coords = data.frame (corJK[, 1], 1))
            # Hotelling correction for n<30
            kk <- ifelse (np <30, min (1, hotelling.transform (kk$corr, kk$dof)), kk$corr)
          }
        }
        cor.matrix <-c (cor.matrix, kk)
      }
    }
  }
  
  # void series with only 1 observation 
  else {
    cor.matrix <- NA
    np <- NA
  }
  # compute mean correlation and sum of obs. pairs
  series.nuts.cor [i] <- mean (cor.matrix, na.rm = T)
  series.nuts.n [i] <- sum (np, na.rm = T)
}

# assemble data frame
intra.nuts.cor <- data.frame ("id" = unique (series.nuts),
                              "cor" = series.nuts.cor,
                              "n" = series.nuts.n)

# plot intra-NUTS correlation
merge.shp.coef <- merge (nutsFasyMap, intra.nuts.cor, by = "id")
final.plot <- merge.shp.coef [order (merge.shp.coef$order), ] 

ggplot () +
  geom_polygon (data = eu, # background admin boundaries
                aes (x = long, y = lat, group = group),
                size = 0.25) + 
  xlim (2500000, 6500000) +
  ylim (1500000, 4500000) +
  geom_polygon (data = final.plot, 
                aes (x = long, y = lat, group = group, fill = cor), 
                color = "black", size = 0.1) +
  scale_fill_gradient2 (guide = "legend", high = "dark red") 
  
# white = Nuts-1 with <2 series or <X years series overlap
# grey = no data

### Nuts-1 chronologies ####

# define function for the mode (with highest values selected if tie)
fmode <- function (x) {
  if (sum (x, na.rm =T)==0) {NA
    } else {
      as.numeric (
        names (which.max (table (x) [order (names (table (x)), decreasing = T)]))
      )
    }
}

# count of sources over Nuts1 by year
df50<-df[df$Yr>=1950,]
df50$Nuts1<-as.character(df50$Nuts1)
no<-data.frame("id"=unique(df50$Nuts1), 'n'=NA)
for (i in 1:nrow(no)) {
  no$n[i]<-length(unique(substr(df50$ID,1,10)[df50$Nuts1==no$id[i]]))
}

# plot beech distribution and number of sources per nuts1

ggplot () +
  geom_polygon (data = eu,
                aes (x = long, y = lat, group = group),
                fill='grey70', color='transparent',size = 0.25) + 
  geom_polygon (data = euforgen, 
                 aes (x = long, y = lat, group = group), 
                 color = NA, fill='grey30',size=0)  +
  xlim (2500000, 5900000) +
  ylim (1500000, 4500000) 

merge.shp.coef <- merge (nutsFasyMap, no, by = "id")
final.plot <- merge.shp.coef [order (merge.shp.coef$order), ] 
final.plot$no.series <- as.factor (final.plot$n)

ggplot () +
  geom_polygon (data = eu, 
                aes (x = long, y = lat, group = group),
                fill='grey70', color='black',size = 0.25) + 
  xlim (2500000, 6500000) +
  ylim (1500000, 4500000) +
  geom_polygon (data = final.plot, 
                aes (x = long, y = lat, group = group, fill = no.series), 
                color = "black", size = 0.1) +
  scale_fill_discrete (guide = "legend", h = c(0, 180) + 15) 

# mode of mast classes over Nuts1 by year
df.nuts.o <- aggregate (df$ORDmast, by=list ("yr" = df$Yr,
                                             "Nuts1" = df$Nuts1), fmode)
colnames (df.nuts.o) <- c ("yr", "nuts", "mi")
df.nuts.o <- df.nuts.o [order (df.nuts.o$yr), ]

# create data frame with masting index by Nuts
series.nuts.o <- reshape (df.nuts.o [, c ("yr", "mi", "nuts")],
                          direction= "wide",
                          idvar = "yr", 
                          timevar = "nuts",
                          v.names = "mi")
names (series.nuts.o) <- gsub ("mi.", "", names (series.nuts.o))
series.nuts.o <- series.nuts.o [,order (names (series.nuts.o))]

# fill empty years
filling.yr <- (min (series.nuts.o$yr):2016) [
  !c (min (series.nuts.o$yr):2016 %in% series.nuts.o$yr)
  ]

filling <- data.frame (matrix (nrow = length (filling.yr), 
                               ncol = ncol (series.nuts.o) - 1), filling.yr)

names (filling) <- names (series.nuts.o)
series.nuts.o <- rbind (series.nuts.o, filling)
series.nuts.o <- series.nuts.o [order (series.nuts.o$yr), ]
series.nuts.oX <- series.nuts.o [, (names (series.nuts.o) != "yr")]
n.nuts <- ncol (series.nuts.oX)
name.nuts <- names (series.nuts.oX)

# Table S1
s1 <- cbind ("yr" = c(1901:2016), series.nuts.oX [series.nuts.o$yr >1900, ])
rownames (s1) <- NULL
#write.csv(s1,"Table_s1.csv")

# define function multiplot
multiplot <- function (..., plotlist = NULL, file, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c (plotlist)
  numPlots <- length (plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null (layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- t (matrix (seq (1, cols * ceiling (numPlots / cols)),
                      ncol = cols, nrow = ceiling (numPlots / cols)))
  }
  if (numPlots == 1) {
    print (plots [[1]])
  } 
  else {
    # Set up the page
    grid.newpage ()
    pushViewport (viewport (layout = grid.layout (nrow (layout), ncol (layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame (which (layout == i, arr.ind = TRUE))
      print (plots [[i]], vp = viewport (layout.pos.row = matchidx$row,
                                         layout.pos.col = matchidx$col))
    }
  }
}

# intersect NUTS-1 and beech distribution
nutsFasyMap <- gIntersection (nuts, euforgen, byid =T, 
                               id = as.character (nuts$NUTS_ID))

names <- data.frame ("NUTS_ID" = row.names (nutsFasyMap))
rownames (names) <- names$NUTS_ID
nutsFasyMap <- SpatialPolygonsDataFrame (nutsFasyMap, data= names)
nutsFasyMap <- fortify (nutsFasyMap, region = "NUTS_ID")

# # color map for years Y1-Y2

# 
#tiff (filename = "Rplot.tiff", width = 1600, height = 1600)

Y1<-1976
Y2<-2014
years<-seq (Y1, Y2)
# cols <- ceiling (sqrt (length (years)))
cols <- 8
plotlist <- list ()

for (i in 1:length (years)) {
  mi <- as.data.frame (t (series.nuts.o [series.nuts.o$yr %in% years [i], ]))
  colnames (mi) <- "index"
  mi$id <- rownames (mi)
  mi <- mi [mi$id != "yr", ]
  merge.shp.coef <- merge (nutsFasyMap, mi, by = "id", all.x = TRUE)
  final.plot <- merge.shp.coef [order (merge.shp.coef$order), ]

  plotlist [[i]] <- ggplot () +
    ggtitle (years [i]) +
    geom_polygon (data = eu,
                  aes (x = long, y = lat, group = group),
                  size = 0, fill='grey90') +
    xlim (2500000, 5900000) +
    ylim (1500000, 4500000) +
    geom_polygon (data = final.plot,
                  aes (x = long, y = lat, group = group, fill = index),
                  color = "transparent", size = 0) +
    scale_fill_gradient (limits = c (1, 5),low = "orange", high = "dark green") +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
}

multiplot (plotlist = plotlist, cols = cols)
#dev.off ()

### temporal autocorrelation (lag 1) ####

# select only 1950-2016 period
initial <- 1950
series.nuts.50X <- series.nuts.oX [series.nuts.o$yr>=initial,]
series.nuts.50 <- series.nuts.o [series.nuts.o$yr>=initial,]

# convert sequences of <X years to NA
series.nuts.50X7 <- series.nuts.50X
for (j in 1:ncol(series.nuts.50X7)) {
  se <- series.nuts.50X7 [, j] 
  seN <- ifelse (is.na (se [1]), NA, 1)
  for (i in 2:length (se)) {
    seN [i] <- ifelse (is.na (se [i]), 0,
                      ifelse (is.na (se [i-1]), 1, seN [i-1] +1))
    }
  seN [is.na (seN)] <- 0 
  
  for (i in 1:X) {
  seN [length(seN)-(i-1)] <- ifelse (seN [length(seN)-0]< (X-i+1), 0, seN [length(seN)-(i-1)])
  }
  
  for (i in length (seN):1) {
    for (k in X:2){
    if (seN[i] ==(k-1)) {seN[i] <- ifelse (seN [i+1] <(k), 0, seN [i])}
    }
  }
  
  seN [seN!=0] <- se [seN!=0]
  seN [seN==0] <- NA
  series.nuts.50X7 [, j] <- seN
}

# prepare data frame
ar <- data.frame ("id"= names (series.nuts.50X), "n" = NA, "ar1" = NA, "p" = NA)

# calculate value for acf coefficient and number of observations
for (i in 1:nrow (ar)) {
  x <- series.nuts.50X7 [, i]
  x1 <- c (NA, x [1:(length (x)-1)])
  df.acf <- na.omit (data.frame ("x" = as.factor (x), x1, "yr"=series.nuts.50$yr))
  ar$n[i] <- nrow (df.acf)
  rm ("x"); rm ("x1")
  if (nrow (df.acf)<2) {ar$ar1[i] <- NA} else 
  if (nrow (df.acf) >=2) {
    model <- clm (x ~ x1, data = df.acf, link = "probit", 
                  control = list (maxIter = 2000))
    ar$ar1 [i] <- round (model$beta, 2)
    sig <- summary(model)$coef ["x1", 4]
    ar$p [i] <- ifelse (is.na (sig), NA, round (sig, 2))
    # create residual chronologies
    # commented out to disable pre-whitening
    # if (!is.na(ar$p [i])) {
    #   if (ar$p [i] <= 0.05) {
    #     df.acf$resid <- as.numeric (as.character (df.acf$x)) - model$fitted.values*5
    #     series.nuts.50X [, i] <- merge (
    #       data.frame ("yr" = series.nuts.50$yr), 
    #       df.acf [, c ("yr", "resid")],
    #       by = "yr", all.x = T) [, "resid"]
    #   }
    # }
  }
}

# check and remove linear trends
ar$slope <- NA
for (i in 1:nrow (ar)) {
  if (ar$n [i] > 1) {
    model <- lm (series.nuts.50X [, i] ~ series.nuts.50$yr)
    ar$slope [i] <- model$coefficients [2]
    if (summary (model) $coefficients [2,4] > 0.05) {ar$slope [i] <- NA}
    if (!is.na (ar$slope [i])) {
      series.nuts.50X [, i] <-  series.nuts.50X [, i] - predict (
        model, newdata = data.frame ("yr" = series.nuts.50$yr))
    }
  }
}

# Table s2
s2 <- ar
#write.csv(s2,"TableS2.csv")

### Mantel test ####

# correlation matrix 
spear.o <- cor (series.nuts.50X, use = "pairwise.complete.obs", method = "spearman")

# label pairs that do not overlap for <X years
cor.controlX <-matrix (NA, n.nuts, n.nuts)
rownames (cor.controlX) <- name.nuts
colnames (cor.controlX) <- name.nuts

for (i in name.nuts){
  for (j in name.nuts){
    corX <- na.omit (merge (series.nuts.50 [, c("yr", i)], 
                            series.nuts.50 [, c("yr", j)], 
                            by="yr"))
    
    if (nrow (corX) < X) {cor.controlX [i, j] = 0}
    if (nrow (corX) >=X) {
      corX <- corX [order (corX$yr), ]
      corX$filter1 <- 0
      corX$filter2 <- 0
      for (k in 1:(max (1, nrow (corX)-X))){
        if (corX$yr [k+X-1] - corX$yr [k] == (X-1)) {
          corX [k, "filter1"] = 1}
      }
      for (k in X:nrow (corX)){
        if (corX$yr [k] - corX$yr [k - (X-1)] == (X-1)) {
          corX [k, "filter2"] = 1}
      }
      corX$filter <- apply (corX [, c ("filter1", "filter2")], 1, max, na.rm=T)
      cor.controlX [i, j] <- sum (corX$filter, na.rm = T)
    }
  }
}

# exclude pairs that do not overlap for at least X years
spear.o [which (cor.controlX < 1)] <- NA

# delete diagonal
for (i in name.nuts) {spear.o [i, i] <- NA}

# simplify similarity matrices by removing NA columns and rows
include <- names (which (apply (spear.o, 2, sum, na.rm = T) != 0))
spear.o <- spear.o [rownames (spear.o) %in% include, colnames (spear.o) %in% include]

# suzuki dissimilarity index
d.o <- as.dist ((-spear.o + 1) / 2 * 100)

# select mapped Nuts1 that overlap for X years
nutsFasySet <- nuts [nuts$NUTS_ID %in% rownames(spear.o), ]
coord <- getSpPPolygonsLabptSlots (nutsFasySet)
coord.labels <- nutsFasySet$NUTS_ID

# distance matrix
dgeo <- dist (coord)

# plot similarity vs distance (figure 2b)
dgeo.m <- as.matrix (dist (coord, diag=T, upper=T))
controprova <- data.frame ("dgeo" = as.vector (dgeo.m),
                           "similarity" = as.vector (spear.o))
controprova <- controprova [order (controprova$dgeo, decreasing = T),]
plot (controprova$dgeo/1000, controprova$similarity, col="dark grey",
                pch=16, cex=.5, ylab = "", xlab = "distance [km]", las=1)
abline (h = 0, col = "grey")
lines(scatter.smooth(controprova$dgeo/1000, controprova$similarity),lwd=3)

# mantel plots and tests
ncf::mantel.test (M1 = as.matrix (dgeo), 
                  M2 = as.matrix (d.o), 
                  resamp = 2000, quiet = T) [1:2] #0.56***

plot.mantel <- mantel.correlog (as.matrix (dgeo), 
                                as.matrix (as.dist (spear.o)), 
                                increment = 400000, 
                                resamp = 2000, 
                                quiet = T)

plot (plot.mantel$mean.of.class/1000, plot.mantel$correlation, 
      xlim = c (0, 2500),
      xlab = "distance [km]",
      ylab = "Spearman's correlation", type = "b", 
      pch = ifelse ((is.na(plot.mantel$p) | plot.mantel$p>0.05), 1, 16), las =1)
abline (h = 0, col = "grey")

# redo with distance along latitude ad longitude only
# by fixing the other coordinate to mean value
coordLat <- coord
coordLat [, 1] <- mean (coordLat [, 1])
dgeoLat <- dist (coordLat)

coordLong <- coord
coordLong [, 2] <- mean (coordLong [, 2])
dgeoLong <- dist (coordLong)

ncf::mantel.test (M1 = as.matrix (dgeoLat), M2 = as.matrix (d.o), resamp = 2000, quiet = T) [1:2] #0.35**
ncf::mantel.test (M1 = as.matrix (dgeoLong),M2 = as.matrix (d.o), resamp = 2000, quiet = T) [1:2] #0.44**

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

### climate correlations ####

# Create lags -1 and -2 year
pre1 <- rbind (clim.na[1:12,], head (pre, (nrow (pre)-12)))
tmx1 <- rbind (clim.na[1:12,], head (tmx, (nrow (tmx)-12)))
tmn1 <- rbind (clim.na[1:12,], head (tmn, (nrow (tmn)-12)))
tmp1 <- rbind (clim.na[1:12,], head (tmp, (nrow (tmp)-12)))
spi31  <- rbind (clim.na[1:12,], head (spi3, (nrow (spi3)-12)))
spei31 <- rbind (clim.na[1:12,], head (spei3, (nrow (spei3)-12)))

pre2 <- rbind (clim.na[1:12,], head (pre1, (nrow (pre1)-12)))
tmx2 <- rbind (clim.na[1:12,], head (tmx1, (nrow (tmx1)-12)))
tmn2 <- rbind (clim.na[1:12,], head (tmn1, (nrow (tmn1)-12)))
tmp2 <- rbind (clim.na[1:12,], head (tmp1, (nrow (tmp1)-12)))
spi32  <- rbind (clim.na[1:12,], head (spi31, (nrow (spi31)-12)))
spei32 <- rbind (clim.na[1:12,], head (spei31, (nrow (spei31)-12)))

pre1$year <- pre2$year <- pre$year
tmp1$year <- tmp2$year <- tmp$year
tmx1$year <- tmx2$year <- tmx$year
tmn1$year <- tmn2$year <- tmn$year
spi31$year <- spi32$year <- spi3$year
spei31$year <- spei32$year <- spei3$year

# change name of months -1 and -2?
pre1$month <- pre$month +12
tmp1$month <- tmp$month +12
tmx1$month <- tmx$month +12
tmn1$month <- tmn$month +12
spi31$month <- spi3$month +12
spei31$month <- spei3$month +12

pre2$month <- pre1$month +12
tmp2$month <- tmp1$month +12
tmx2$month <- tmx1$month +12
tmn2$month <- tmn1$month +12
spi32$month <- spi31$month +12
spei32$month <- spei31$month +12

# combine the datasets from lag 0 to -2 
pre.all <- rbind (pre, pre1, pre2)
tmn.all <- rbind (tmn, tmn1, tmn2)
tmx.all <- rbind (tmx, tmx1, tmx2)
tmp.all <- rbind (tmp, tmp1, tmp2)
spi3.all <- rbind (spi3, spi31, spi32)
spei3.all <- rbind (spei3, spei31, spei32)

# subset the climate data to include only NUTS with mast data
pre.all <- cbind (pre.all [, c ("year", "month")], pre.all [, names (pre.all) %in% nutsDef])
tmn.all <- cbind (tmn.all [, c ("year", "month")], tmn.all [, names (tmn.all) %in% nutsDef])
tmp.all <- cbind (tmp.all [, c ("year", "month")], tmp.all [, names (tmp.all) %in% nutsDef])
tmx.all <- cbind (tmx.all [, c ("year", "month")], tmx.all [, names (tmx.all) %in% nutsDef])
spi3.all <- cbind (spi3.all [, c ("year", "month")], spi3.all [, names (pre.all) %in% nutsDef])
spei3.all <- cbind (spei3.all [, c ("year", "month")], spei3.all [, names (spei3.all) %in% nutsDef])

# subset the (detrended) masting data
mast <- series.nuts.50X [series.nuts.50$yr >= initial, ]

# count non-NA observations in the study period
nuts.obs <- apply (mast, 2, function (x) {length (na.omit (x))})

# exlude NUTS with fewer than X years of observations
mast <- mast [, which (nuts.obs >= X)]

# Note that depending on the period of analysis, there may be issues with some other NUTS that have >7 obs, but few non-mast years
# E.g. DE8 has greater than 8 individual years, but they are all CAT 5 years
mast$DE8<- NULL

# reorder masting data by NUTS1 with year at the beginning
mast <- data.frame ("yr" = c (initial:2016), mast)

# function for calculating the spearman rank coeff
myfunction <- function(var1, var2, thresh=X){
  AA<- cor.test(var1, var2, use="complete.obs", method="spearman")$estimate
  A<- ifelse(length(which(var1 != "NA")) >=thresh, AA, NA)
  return(A)
}

# function for calculating the spearman rank p-value
myfunction2 <- function(var1, var2, thresh=X){
  AA<- cor.test(var1, var2, use="complete.obs", method="spearman")$p.value
  A<- ifelse(length(which(var1 != "NA")) >=thresh, AA, NA)
  return(A)
}

# plotting functions
# plots
col20<- rgb(0, 20, 200, 255, maxColor=255) # blue end
col19<- rgb(0, 84, 255, 255, maxColor=255)
col17<- rgb(0, 168, 255, 255, maxColor=255) 
col13<- rgb(230, 250, 255, 255, maxColor=255)
col11<- rgb(255, 255, 255, 255, maxColor=255) # white
col9<- rgb(255, 255, 220, 255, maxColor=255)
col4<- rgb(255, 150, 0, 255, maxColor=255)
col2<- rgb(255, 60, 0, 255, maxColor=255) 
col1<- rgb(200, 20, 0, 255, maxColor=255) # red end

# function for creating a colour bar
color.bar <- function (lut, min, max=-min, nticks=18, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, lwd=0, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# Some plotting parameters
m<- ncol(mast)   #number of NUTs regions
n<- 36    # number of months
x <- seq(m)
y <- seq(n)

ylabs<- c("DEC", "NOV", "OCT", "SEP", "AUG", "JUL", "JUN", "MAY", "APR", "MAR", "FEB", "JAN", 
          "DEC.1", "NOV.1", "OCT.1", "SEP.1", "AUG.1", "JUL.1", "JUN.1", "MAY.1", "APR.1", "MAR.1", "FEB.1", "JAN.1",
          "DEC.2", "NOV.2", "OCT.2", "SEP.2", "AUG.2", "JUL.2", "JUN.2", "MAY.2", "APR.2", "MAR.2", "FEB.2", "JAN.2")

# create colour scale
rbPal <- colorRampPalette (c (col20, col19, col17, col13,col11, col9, col4, col2, col1))
Col<- rbPal(100)

# PLOT
nNUTS <- ncol (mast) -1 #number of nuts regions
nMonths <- n #number of months to be investigated

spearman.detrend <- list ()
p.value.detrend <- list ()

clim <- list (pre.all, tmn.all, tmx.all, tmp.all, spi3.all, spei3.all)
clim.names<- c("pre.all", "tmn.all", "tmx.all", "tmp.all", "spi3.all", "spei3.all")


# loop function to calculate the spearman ranks for each NUTS and each month
# this is a bit slow...

for (k in 1:length (clim)) {

  spearman.detrend [[k]] <- data.frame (matrix (NA, nrow = nMonths, ncol = nNUTS))
  p.value.detrend [[k]] <- data.frame (matrix (NA, nrow = nMonths, ncol = nNUTS))
  for (i in 1:nMonths)  # Month
  {
    for(j in 1:nNUTS) # NUTS
    {
      clim[[k]] <- clim[[k]][clim[[k]]$year>=initial,] 
      a <- (subset (clim[[k]], clim[[k]][,"month"] == i)[,j+2])
      b <- (1:length (subset (clim[[k]], clim[[k]][,"month"] == i)[,j+2]))
      A <- lm (a~b)
      clim_resid <- data.frame (rep (NA, length (mast [,j])))
      names (clim_resid) <- "residual"
      clim_resid [names (A$residuals), "residual"] <- A$residuals
      spearman.detrend [[k]][i,j] <- myfunction (mast[,j+1], clim_resid$residual)
      p.value.detrend [[k]] [i,j] <- myfunction2 (mast[,j+1], clim_resid$residual)
    }
  }

  # then seasonal climate
  
  spearman.seasonal <- list ()
  p.value.seasonal <- list ()
  
  spearman.seasonal [[k]] <- data.frame (matrix (NA, nrow = 3, ncol = nNUTS))
  p.value.seasonal [[k]] <- data.frame (matrix (NA, nrow = 3, ncol = nNUTS))
  for (i in 1:3)  
  {
    for(j in 1:nNUTS) # NUTS
    {
      clim[[k]] <- clim[[k]][clim[[k]]$year>=initial,] 
      
      jj.1 <- (subset (clim[[k]], clim[[k]][,"month"] == 18)[,j+2]) + (subset (clim[[k]], clim[[k]][,"month"] == 19)[,j+2])
      jj.2 <- (subset (clim[[k]], clim[[k]][,"month"] == 30)[,j+2]) + (subset (clim[[k]], clim[[k]][,"month"] == 31)[,j+2])
      dT<- jj.1-jj.2
      
      jj.1.yr<- (1:length (subset (clim[[k]], clim[[k]][,"month"] == 18)[,j+2]))
      jj.2.yr<- (1:length (subset (clim[[k]], clim[[k]][,"month"] == 30)[,j+2]))
      dT.yr<- (1:length (dT))
      
      MOD.jj.1 <- lm (jj.1~jj.1.yr)
      jj.1_resid <- data.frame (rep (NA, length (mast [,j])))
      names (jj.1_resid) <- "residual"
      jj.1_resid[names (MOD.jj.1$residuals), "residual"] <- MOD.jj.1$residuals
      
      MOD.jj.2 <- lm (jj.2~jj.2.yr)
      jj.2_resid <- data.frame (rep (NA, length (mast [,j])))
      names (jj.2_resid) <- "residual"
      jj.2_resid[names (MOD.jj.2$residuals), "residual"] <- MOD.jj.2$residuals
      
      MOD.dT <- lm (dT~dT.yr)
      dT_resid <- data.frame (rep (NA, length (mast [,j])))
      names (dT_resid) <- "residual"
      dT_resid[names (MOD.dT$residuals), "residual"] <- MOD.dT$residuals
      
      spearman.seasonal [[k]][3,j] <- myfunction (mast[,j+1], jj.2_resid$residual)
      p.value.seasonal[[k]] [3,j] <- myfunction2 (mast[,j+1], jj.2_resid$residual)
      
      spearman.seasonal [[k]][2,j] <- myfunction (mast[,j+1], jj.1_resid$residual)
      p.value.seasonal[[k]] [2,j] <- myfunction2 (mast[,j+1], jj.1_resid$residual)
      
      spearman.seasonal [[k]][1,j] <- myfunction (mast[,j+1], dT_resid$residual)
      p.value.seasonal[[k]] [1,j] <- myfunction2 (mast[,j+1], dT_resid$residual)
      
    }
  }
  
  # then select the values for plotting
  data <- spearman.detrend [[k]]
  p.value <- p.value.detrend [[k]]
  data.s <- spearman.seasonal [[k]]
  p.value.s <- p.value.seasonal [[k]]
  
  colnames (data)<- c (as.character (colnames (mast[2:ncol(mast)])))
  colnames (p.value)<- c (as.character (colnames (mast[2:ncol(mast)])))
  colnames (data.s)<- c (as.character (colnames (mast[2:ncol(mast)])))
  colnames (p.value.s)<- c (as.character (colnames (mast[2:ncol(mast)])))
  
  # re-order according to cluster group
  # extract NC from Cluster 1
  
  cluster.1<- mycl.o$id [which (mycl.o$cluster == 1)]
  names.cluster.1<- match (cluster.1, colnames (data))
  data.cluster.1 <- data [,names.cluster.1] 
  p.value.cluster.1 <- p.value [,names.cluster.1]
  data.cluster.1.s <- data.s [,names.cluster.1] 
  p.value.cluster.1.s <- p.value.s [,names.cluster.1]
  
  # extract NC from Cluster 2
  cluster.2<- mycl.o$id [which (mycl.o$cluster == 2)]
  names.cluster.2<- match (cluster.2, colnames (data))
  data.cluster.2 <- data [,names.cluster.2] 
  p.value.cluster.2 <- p.value [,names.cluster.2]
  data.cluster.2.s <- data.s [,names.cluster.2] 
  p.value.cluster.2.s <- p.value.s [,names.cluster.2]
  
  # extract NC from Cluster 3
  cluster.3<- mycl.o$id [which(mycl.o$cluster == 3)]
  names.cluster.3<- match (cluster.3, colnames (data))
  data.cluster.3 <- data [,names.cluster.3] 
  p.value.cluster.3 <- p.value [,names.cluster.3]
  data.cluster.3.s <- data.s[,names.cluster.3] 
  p.value.cluster.3.s <- p.value.s [,names.cluster.3]
  
  # some regions are not in the cluster plot, but are included in the correlation analysis
  data.cluster.none <- subset (data, select = -c(names.cluster.1, names.cluster.2, names.cluster.3))
  p.value.cluster.none <- subset (p.value, select = -c(names.cluster.1, names.cluster.2, names.cluster.3))
  data.cluster.none.s <- subset (data.s, select = -c(names.cluster.1, names.cluster.2, names.cluster.3))
  p.value.cluster.none.s <- subset (p.value.s, select = -c(names.cluster.1, names.cluster.2, names.cluster.3))
  
  # data for plotting
  # recombine the data according to clusters
  data.plot<- cbind (data.cluster.1, data.cluster.2, data.cluster.3, data.cluster.none)
  p.value.plot<- cbind (p.value.cluster.1, p.value.cluster.2, p.value.cluster.3, p.value.cluster.none)
  
  data.plot.s<- cbind (data.cluster.1.s, data.cluster.2.s, data.cluster.3.s, data.cluster.none.s)
  p.value.plot.s<- cbind (p.value.cluster.1.s, p.value.cluster.2.s, p.value.cluster.3.s, p.value.cluster.none.s)
  
  xlabs<- c(as.character(colnames(data.plot)))
  xlabs2<- apply(mast, 2, function(x) length(which(!is.na(x))))
  xlabs2<- data.frame(xlabs2[2:length(xlabs2)])
  xlabs2<- xlabs2[xlabs,]
  
  #------#
  
  png(paste("FINAL_MATRIX-PLOT_", clim.names[k], ".png", sep=""), width=21, height=12, units="cm", res=800)
  
  par(mar=c(3,3,0.6,3)+0.1, mgp=c(1.1,0.2,0), tck=0.02, lwd=1, cex=1, cex.lab=0.5, cex.axis=0.5,  pch=16)
  par(xaxs="i"); par(yaxs="i")
  par(fig=c(0,0.96,0.055,1))
  
  # plot the correlation values according to the colour scale
  # plotted according to order JAN-2 --> DEC0
  z <-data.matrix(rbind(data.plot[12:1,], data.plot[24:13,], data.plot[36:25,])) #rownames.force =NA)
  image(x,y,t(z), col= Col, axes=FALSE, ylab="", xlab="", zlim=c(-1, 1))
  a<-seq(1:n)
  dataA<- matrix(c(rep(a,nNUTS)), n, nNUTS, byrow=F)
  a<-seq(1:nNUTS)
  dataB<- matrix(c(rep(a,n)), n, nNUTS, byrow=T)
  
  # plot p-values
  z.p <-data.matrix(rbind(p.value.plot[12:1,], p.value.plot[24:13,], p.value.plot[36:25,])) #rownames.force =NA)
  points(dataB+0.5, dataA, col="black", pch=ifelse(z.p < (0.01/39)*5, "*", NA_integer_), cex=1, font=2) #rgb=c(5,5,5,100, maxColorValue=255)
  points(dataB+0.5, dataA, col="black", pch=ifelse(z.p < 0.01*5, 16, NA_integer_), cex=0.2) #rgb=c(5,5,5,100, maxColorValue=255)
  
  #mtext("NUTS REGION",1, line=1, cex=0.7)
  mtext("MONTH",2, line=1.5, cex=0.7, las=0)
  axis(2, at=c(1:n), tcl=0, labels=ylabs, las=1)
  
  COLS.NUTS<- c(rep("coral1", length(cluster.1)),rep("green3", length(cluster.2)),rep("dodgerblue", length(cluster.3)),rep("black", ncol(data.cluster.none)))
  #axis(3, at=c(1:m), tcl=0, labels=, las=0, cex=0.5, col.axis=COLS.NUTS)
  text(1:nNUTS+0.5, rep(n+1, nNUTS), xlabs2, xpd=TRUE, cex=0.5, col=COLS.NUTS)
  
  grid(nx = nNUTS, ny=n, col="white", lty="solid", lwd=2)
  box(col="white", lwd=2)
  box()
  
  # then add a column for the mean CC
  par(mar=c(3,1,.6,1)+0.1, mgp=c(1.1,0.2,0), tck=0.02, lwd=1, cex=1, cex.lab=0.5, cex.axis=0.5,  pch=16)
  par(xaxs="i")
  par(yaxs="i")
  par(new=TRUE)
  par(fig=c(0.868,0.941,0.055,1))
  
  zz <-data.matrix(rowMeans(data, na.rm=TRUE)*2) # curremntly multiplied by 2 for plotting
  z2 <-data.matrix(c(zz[12:1,], zz[24:13,], zz[36:25,]))
  
  image(1,y,t(z2), col=Col, axes=FALSE, ylab="", xlab="", zlim=c(-1, 1))
  #mtext("MEAN*",1, line=0.1, cex=0.5, las=2)
  #mtext("*X2",1, line=1.1, cex=0.45, las=0) # needs adjusting to reflect above
  
  grid(nx = 1, ny=n, col="white", lty="solid", lwd=2)
  box(col="white", lwd=6)
  box()
  
  # then add a row for the seasonal CC and deltaT
  par(new=TRUE) 
  par(mar=c(2.5,3,0.4,3)+0.1, mgp=c(1.1,0.2,0), tck=0.02, lwd=1, cex=1, cex.lab=0.5, cex.axis=0.5,  pch=16)
  par(xaxs="i")
  par(yaxs="i")
  par(fig=c(0,0.96, 0, 0.195))
  
  z <-data.matrix(data.plot.s[1:3,]) #rownames.force =NA)
  image(x,1:3,t(z), col= Col, axes=FALSE, ylab="", xlab="", zlim=c(-1, 1))
  
  a<-seq(1:3)
  dataA<- matrix(c(rep(a,nNUTS)), 3, nNUTS, byrow=F)
  a<-seq(1:nNUTS)
  dataB<- matrix(c(rep(a,n)), 3, nNUTS, byrow=T)
  
  # plot p-values
  z.p <-data.matrix(p.value.plot.s[1:3,]) #rownames.force =NA)
  points(dataB+0.5, dataA, col="black", pch=ifelse(z.p < (0.01*5)/39, "*", NA_integer_), cex=1, font=2) #rgb=c(5,5,5,100, maxColorValue=255)
  points(dataB+0.5, dataA, col="black", pch=ifelse(z.p < 0.01*5, 16, NA_integer_), cex=0.2) #rgb=c(5,5,5,100, maxColorValue=255)
  
  mtext("NUTS REGION",1, line=1, cex=0.7)
  #mtext("MONTH",2, line=1.5, cex=0.7, las=0)
  axis(2, at=c(1:3), tcl=0, labels=c("dT", "JJ.1", "JJ.2"), las=1)
  
  COLS.NUTS<- c(rep("coral1", length(cluster.1)),rep("green3", length(cluster.2)),rep("dodgerblue", length(cluster.3)),rep("black", ncol(data.cluster.none)))
  axis(1, at=c(1:m), tcl=0, labels=FALSE, las=2)
  text((1:nNUTS)+0.5, rep(-0.7, nNUTS), xlabs, xpd=TRUE, cex=0.5, srt=90, col=COLS.NUTS)
  
  grid(nx = nNUTS, ny=3, col="white", lty="solid", lwd=2)
  box(col="white", lwd=2)
  box()
  
  # then add a column for the mean seasonal CC
  par(mar=c(2.5,1,0.4,1)+0.1, mgp=c(1.1,0.2,0), tck=0.02, lwd=1, cex=1, cex.lab=0.5, cex.axis=0.5,  pch=16)
  par(xaxs="i")
  par(yaxs="i")
  par(new=TRUE)
  par(fig=c(0.868,0.941,0,0.195))
  
  zz <-data.matrix(rowMeans(data.s, na.rm=TRUE)*2) # curremntly multiplied by 2 for plotting
  image(1,1:3,t(zz), col=Col, axes=FALSE, ylab="", xlab="", zlim=c(-1, 1))
  mtext("MEAN*",1, line=0.1, cex=0.5, las=2)
  mtext("*X2",1, line=1.1, cex=0.45, las=0) # needs adjusting to reflect above
  
  grid(nx = 1, ny=3, col="white", lty="solid", lwd=2)
  box(col="white", lwd=2)
  abline(v=c(0.56), lwd=8, col="white")
  abline(v=c(1.44), lwd=8, col="white")
  
  box()
  
  # add the colour scale
  
  par(mar=c(0.5,0.5,.4,.5)+0.1, mgp=c(0,0.2,0), tck=0.08, lwd=1, cex=1, cex.lab=0.5, cex.axis=0.5,  pch=16)
  par(xaxs="i")
  par(yaxs="i")
  
  par(new=TRUE)
  par(fig=c(0.94,0.99,0.2,0.8))
  
  max(z2);min(z2)
  
  color.bar(rbPal(100), max=1, min=-1, nticks=5)
  axis(2, at=c(-1.0, -0.5, 0.0, 0.5, 1.0),las=1, labels=FALSE)
  axis(4, at=c(-1, -0.45, 0.00, 0.45, 1),las=1, labels=FALSE, tck=0)
  axis(1,las=1, labels=FALSE, tck=0)
  axis(3,las=1, labels=FALSE, tck=0)
  
  dev.off()
}

# end of correlation plotting code

# correlation between spearman and latitude
rho.max <- data.frame (spearman.detrend [[3]])
rho.pre <- data.frame (spearman.detrend [[1]])
p.max <- data.frame (p.value.detrend  [[3]])

colnames (rho.max) <- names (mast) [-1]
colnames (rho.pre) <- names (mast) [-1]
rownames (rho.max) <- c (ylabs[12:1], ylabs[24:13], ylabs[36:25])
rownames (rho.pre) <- c (ylabs[12:1], ylabs[24:13], ylabs[36:25])

lat <- coordDef [which (nutsDef %in% colnames (rho.max)),2] 
lon <- coordDef [which (nutsDef %in% colnames (rho.max)),2] 

lat.max <- c(); lat.pre <- c()
lon.max <- c(); lon.pre <- c()

for (i in 1:nrow (rho.max)) {
  lat.max [i] <- cor (unlist (rho.max [i,]), lat)
  lat.pre [i] <- cor (unlist (rho.pre [i,]), lat)
  lon.max [i] <- cor (unlist (rho.max [i,]), lon)
  lon.pre [i] <- cor (unlist (rho.pre [i,]), lon)
}

lat.all <- data.frame ("max"=lat.max, "pre"=lat.pre)
lon.all <- data.frame ("max"=lon.max, "pre"=lon.pre)

rownames (lat.all) <- rownames (lon.all) <- rownames (rho.pre)

summerCor <- rho.max[c("JUN.1","JUL.1","AUG.1","JUN.2","JUL.2","AUG.2"),]
summerP <- p.max [rownames(rho.max) %in% c ("JUN.1","JUL.1","AUG.1","JUN.2","JUL.2","AUG.2"),]
summerP [summerP<=0.05]<-1
summerP [summerP>0.05 & summerP<1]<-16
summerP <- as.numeric(summerP)
rownames(summerCor) <- c("JUN-1","JUL-1","AUG-1","JUN-2","JUL-2","AUG-2")

par(mfrow=c(3,3))
for (i in 1:nrow(summerCor)){
  plot (lat, unlist (summerCor [i,]), main=rownames(summerCor)[i], 
                  ylab="R masting vs. Tmax", xlab = "latitude", pch=16, 
                  ylim=c(-0.8,0.8), col=as.numeric(summerP[,i]))
  linear<-predict(lm (unlist (summerCor [i,])~lat),interval="prediction")
  lines(lat,data.frame(linear)$fit)
  lines(lat,data.frame(linear)$lwr,lty=3)
  lines(lat,data.frame(linear)$upr,lty=3)
}

# par(mfrow=c(2,2))
# scatter.smooth (lat, unlist (summerCor [1,])-unlist (summerCor [2,]), type="p", main="June-July(-1)", ylab="Delta R masting vs. Tmax", pch=16)
# scatter.smooth (lat, unlist (summerCor [2,])-unlist (summerCor [3,]), type="p", main="July-August(-1)", ylab=" ", pch=16)
# scatter.smooth (lat, unlist (summerCor [4,])-unlist (summerCor [5,]), type="p", main="June-July(-2)", ylab="Delta R masting vs. Tmax", pch=16)
# scatter.smooth (lat, unlist (summerCor [5,])-unlist (summerCor [6,]), type="p", main="July-August(-2)", ylab=" ", pch=16)
# 
# summary(lm (unlist (summerCor [1,])-unlist (summerCor [2,])~lat ))
# summary (lm(unlist (summerCor [2,])-unlist (summerCor [3,])~lat ))
# 
# par(mfcol=c(1,2))
summerMaxCor1<-as.character(apply(summerCor[1:3,],2,which.is.max))
summerMaxCor2<-as.character(apply(summerCor[4:6,],2,which.is.max))

summerMaxCor1[summerMaxCor1==1] <- "Jun"
summerMaxCor1[summerMaxCor1==2] <- "Jul"
summerMaxCor1[summerMaxCor1==3] <- "Aug"
summerMaxCor2[summerMaxCor2==1] <- "Jun"
summerMaxCor2[summerMaxCor2==2] <- "Jul"
summerMaxCor2[summerMaxCor2==3] <- "Aug"

boxplot(lat~summerMaxCor1, ylab="latitude",xlab="month with max R vs. Tmax",varwidth=T)
boxplot(lat~summerMaxCor2, ylab=" ",xlab="month with max R vs. Tmax",varwidth=T)

par(mfcol=c(1,1))

### Ordinal logistic regression model ####
# NOTES THAT THE FULL CLIMATE DATASET IS REQUIRED (1901-2016)
# IF THE CODE ABOVE WAS RUN FOR THE PERIOD 1950-2016 THEN IT WILL NEED RERUNNING
# WITHOUT TRUNCATING THE CLIMATE AND MAST DATA TO 1950-2016

# run when testing

# mast<- read.table("mast.prelim.delete.txt", header=TRUE)
# pre<- read.csv("pre.2017.csv", header=TRUE)
# tmx<- read.csv("tmx.2017.csv", header=TRUE)

sites <- c ("DE1", "DE2", "DE9", "DEF", "DK0", "NL1", "SE2", "UKJ")
pre.8 <- pre[,c(sites, "anni", "mesi")]
tmx.8 <- tmx[,c(sites, "anni", "mesi")]
mast.8<- mast[,c(sites, "yr")]

clim.na<- as.data.frame (matrix (rep (NA, 12*length (tmx.8)), ncol = length (tmx.8)))
names (clim.na) <- names (tmx.8)

pre1 <- rbind (clim.na, head (pre.8, (nrow (pre.8)-12)))
tmx1 <- rbind (clim.na, head (tmx.8, (nrow (tmx.8)-12)))

pre2 <- rbind (clim.na, clim.na, head (pre.8, (nrow (pre.8)-24)))
tmx2 <- rbind (clim.na, clim.na, head (tmx.8, (nrow (tmx.8)-24)))

# change name of mesis -1 and -2?
pre1$mesi <- pre1$mesi +12
tmx1$mesi <- tmx1$mesi +12

pre2$mesi <- pre2$mesi +24
tmx2$mesi <- tmx2$mesi +24


# change the annis
pre1$anni<- as.numeric(pre1$anni) + 1
tmx1$anni<- as.numeric(tmx1$anni) + 1
pre2$anni<- as.numeric(pre2$anni) + 2
tmx2$anni<- as.numeric(tmx2$anni) + 2

# combine the datasets from lag 0 to -2 
pre.all <- rbind (pre.8, pre1, pre2)
tmx.all <- rbind (tmx.8, tmx1, tmx2)

# remove annis 1901 and 1902 as do not have complete data (no -1 and -2 climate data)
pre.all<- subset(pre.all, pre.all$anni >1902)
tmx.all<- subset(tmx.all, tmx.all$anni >1902)

# creating new datasets
initial<- 1950
initial2 <- 1901
mast <- series.nuts.oX [series.nuts.o$yr >= initial2, ]
mast <- data.frame ("yr" = c (initial2:2016), mast)
# also remove 1901 and 1902 from mast data
mast<- subset(mast, mast$yr > 1902)

sites <- c ("DE1", "DE2", "DE9", "DEF", "DK0", "NL1", "SE2", "UKJ")
clim <- list (pre.all, tmx.all)

# check and remove linear trends
slope <- NA
for (i in 1:ncol (series.nuts.oX)) {
  if (length (na.omit (mast [, i+1])) <3) {slope[i] <- NA}
  else {
    model <- lm (mast [, i+1] ~ mast$yr)
    slope [i] <- model$coefficients [2]
    if (summary (model) $coefficients [2,4] > 0.05) {slope [i] <- NA}
  }
} 
slope <- cbind (slope, "id" = names (series.nuts.oX))
slope # no slope for 8 longest NC

clim.ord <- data.frame (clim[[1]] [clim[[1]]$anni>=initial2, c("anni", "mesi", sites)])
clim2.ord <- data.frame (clim[[2]] [clim[[1]]$anni>=initial2, c("anni", "mesi", sites)])
mast.ord <- data.frame (mast [,c("yr", sites)])


JUL.1<- 19 # mesis used for modelling
JUL.2<- 31
JUN.1<- 18 # mesis used for modelling
JUN.2<- 30

# Building models for each NUTS
coefs.all_models <- matrix(nrow=13, ncol=8, byrow=FALSE)
stats.all_models <- matrix(nrow=14, ncol=8, byrow=FALSE)
validation.all_models<- matrix(nrow=1, ncol=8, byrow=FALSE)
k_fold.val<- matrix(nrow=1, ncol=8, byrow=FALSE)
# NOTE THIS A VERY LONG FUNCTION!

for(i in 1:8) { # sites
  site<- sites[i]
  # first, create the detrended climate indices (8 in total, 4 TMAX and 4 PRE)
  # precip
  mesi<- 7+12 #JUL.1
  a<- scale((clim[[1]][clim[[1]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[1]][clim[[1]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  PREC.JUL.1<- clim_resid[,1]
  
  mesi<- 6+12 #JUN.1
  a<- scale((clim[[1]][clim[[1]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[1]][clim[[1]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  PREC.JUN.1<- clim_resid[,1]
  
  mesi<- 7+24 #JUL.2
  a<- scale((clim[[1]][clim[[1]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[1]][clim[[1]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  PREC.JUL.2<- clim_resid[,1]
  
  mesi<- 6+24 #JUN.2
  a<- scale((clim[[1]][clim[[1]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[1]][clim[[1]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  PREC.JUN.2<- clim_resid[,1]
  
  #temp
  mesi<- 7+12 #JUL.1
  a<- scale((clim[[2]][clim[[2]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[2]][clim[[2]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  TEMP.JUL.1<- clim_resid[,1]
  
  mesi<- 6+12 #JUN.1
  a<- scale((clim[[2]][clim[[2]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[2]][clim[[2]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  TEMP.JUN.1<- clim_resid[,1]
  
  mesi<- 7+24 #JUL.2
  a<- scale((clim[[2]][clim[[2]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[2]][clim[[2]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  TEMP.JUL.2<- clim_resid[,1]
  
  mesi<- 6+24 #JUN.2
  a<- scale((clim[[2]][clim[[2]][,10] == mesi,site]))[,1]
  b<- (1:length(clim[[2]][clim[[2]][,10] == mesi,site]))
  A<- lm(a~b)
  clim_resid<- data.frame(rep(NA, nrow(mast)))
  names(clim_resid)<- "residual"
  clim_resid[names(A$residuals), "residual"]<- A$residuals
  TEMP.JUN.2<- clim_resid[,1]
  
  # calculate previous anni mast value (for NC-1)
  MAST.L<- scale(c(NA, mast[,site][1:nrow(mast)-1]))[,1]
  
  # create data for model
  data<- data.frame (cbind (mast.ord [,site], MAST.L, TEMP.JUL.1, TEMP.JUN.1, TEMP.JUL.2, TEMP.JUN.2, PREC.JUL.1, PREC.JUN.1, PREC.JUL.2, PREC.JUN.2))
  data$anni<- c(1903:2016)
  
  # split into calibration and validation datasets
  data.cc<- subset(data, data$anni>=initial)
  data.vv<- subset(data, data$anni<initial)
  # remove NA values
  data.cc<- data.frame(data.cc[complete.cases(data.cc),])
  data.vv<- data.frame(data.vv[complete.cases(data.vv),])
  
  # TEMP.JUL.1 and PREC.JUL.1 frequently colinear, so remove PREC.JUL.1
  # plot(TEMP.JUL.1~PREC.JUL.1)
  
  # now fit automatically using stepAIC to find the optimal model 
  # saturated model
  fit1 <- lrm(data.cc[,1] ~  MAST.L + 
                TEMP.JUL.1 + 
                TEMP.JUL.2 + 
                TEMP.JUN.1 +
                TEMP.JUN.2 + 
                PREC.JUL.2 + 
                PREC.JUN.1 + 
                PREC.JUN.2, 
              data=data.cc)
  
  # fastbw() the fitting function from the rms package
  step<- fastbw(fit1, rule="aic", type="individual")
  
  # keep the selected variables
  select.var<- data.cc[,step$names.kept]
  
  # fit the model with with the selected variables
  optimal<- lrm.fit(select.var, data.cc[, 1])
  optimal.mod<- lrm(data.cc[,1]~ ., data = select.var, x=TRUE, y=TRUE)
  
  # create a dataframe with model coefs, stats and validation
  coefs.master<- data.frame(names=c("y>=2", "y>=3", "y>=4", "y>=5", "MAST.L", "TEMP.JUL.1", "TEMP.JUN.1" ,"TEMP.JUL.2", "TEMP.JUN.2","PREC.JUL.1", "PREC.JUN.1" ,"PREC.JUL.2", "PREC.JUN.2"), rep(NA, 13))
  coefs.optimal<- data.frame(optimal$coefficients)
  coefs.optimal$names<- as.character(rownames(coefs.optimal))
  coefs<- (merge(coefs.master,coefs.optimal,by="names", all.x=TRUE))
  coefs<- coefs[order(coefs$names),][,"optimal.coefficients"]
  
  stats.optimal<- data.frame(optimal$stats)[,1]
  
  coefs.all_models[,i]<- coefs
  stats.all_models[,i]<- stats.optimal
  
  # validation
  model.predict<- predict(optimal.mod, type=c("mean"), newdata=data.vv, se.fit=FALSE, codes=FALSE)
  observed<- data.vv[,"V1"]
  validation.all_models[,i]<- summary(lm(model.predict~ observed))$r.squared
  
  # LOOCV
  set.seed(15)
  k_fold <- validate(optimal.mod, method="boot", B=1000,
                     bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0,
                     force=NULL, estimates=TRUE, pr=FALSE, group=data.cc[,1])
  k_fold.val[i] <- k_fold[2,5]
  
}

# give objects column names and row names
colnames(coefs.all_models)<- sites
colnames(stats.all_models)<- sites
colnames(validation.all_models)<- sites
colnames(k_fold.val)<- sites

rownames(coefs.all_models)<- c("MAST.L", "PREC.JUL.1", "PREC.JUL.2", "PREC.JUN.1", "PREC.JUN.2", "TEMP.JUL.1", "TEMP.JUL.2", "TEMP.JUN.1", "TEMP.JUN.2", "y>=2", "y>=3", "y>=4", "y>=5")
rownames(stats.all_models)<- names(optimal$stats)

# write datasets if required
# write.table(coefs.all_models, "Model_coefs.txt")
# write.table(stats.all_models, "Model_stats.txt")
# write.table(validation.all_models, "Model_validation(Rho).txt")

# plot comparing the different NUTS models

#tiff("Figure_5_revised.tif", width=8, height=12, units="cm", res=800)

rbPal <- colorRampPalette(c(col20, col17, "white", col4,col1))

par(mar=c(1,2.7,.1,.5)+0.1, mgp=c(1.5,0.2,0), tck=0.08, lwd=1, cex=1, cex.lab=0.6, cex.axis=0.6,  pch=16, tck=0)
par(fig=c(0,1,0.6,1))

plot(1:8-0.2, stats.all_models["R2",], axes=F, pch=15, type="h", lwd=1.5, ylab= "", 
     xlab="", col=ifelse(stats.all_models["P",]> 0.05, "grey80", "black"), ylim=c(0,0.87), xlim=c(0.85, (8+0.2)))
points(1:8, k_fold.val[1,], type="h", lwd="1.5", cex=0.9, col=ifelse(stats.all_models["P",]> 0.05, "grey90", "grey60"))
points(1:8+0.2, validation.all_models[1,], type="h", lwd="1.5", cex=0.9, col=ifelse(stats.all_models["P",]> 0.05, "grey90", "grey75"))

points(1:8-0.2, stats.all_models["R2",], pch=16, cex=0.8, col=ifelse(stats.all_models["P",]> 0.05, "grey90", "black"))
points(1:8, k_fold.val[1,], pch=16, cex=0.8, col=ifelse(stats.all_models["P",]> 0.05, "grey90", "grey60"))
points(1:8+0.2, validation.all_models[1,], pch=16, cex=0.8, col=ifelse(stats.all_models["P",]> 0.05, "grey90", "grey75"))

axis(2, at=c(0:8/10), labels=c(0:8/10), las=1, tck=0.01)
mtext("Model statistics", 2, line=1, cex=0.8)
ylab<- colnames(coefs.all_models)
axis(1, at=c(1:8), labels=FALSE, las=2)
text(1:8, rep(-0.11, 8), ylab, xpd=TRUE, srt=90, cex=0.6, col=c("black"))
mtext(text="A)", side=2,  padj=-6.5, las=1, line=1.5, cex=0.9)

legend("top", cex=0.5, xjust=0.5, c("Rsq_1951-2015", "Rsq_bootstrap", "Rsq_1901-1950"), col=c("black", "grey60", "grey70"),  pch=16, ncol=3, bty="n", x.intersp = 0.4)
box()

par(new=TRUE)
par(fig=c(0,1,0.51,0.59))
image(x=1:8, y=1, matrix(rep(0,8)), zlim=c(-9, 9), col="white", axes=F, ylab="", xlab="")
axis(2, at=c(1), labels="OBS", las=1)
text((1:8), (rep(1, 9)), c(stats.all_models["Obs",]), cex=0.6)
box()
abline(v=(1:13)+0.5)

par(new=TRUE)
par(fig=c(0,1,0.08,0.55))
coefs.matrix<- matrix(nrow=9, ncol=8, coefs.all_models[c(1,6,8,7,9,2,4,3,5),])
Col<- rbPal(100)
image(x=1:8, y=1:9, t(coefs.matrix), zlim=c(-2.5, 2.5), col=Col, ylim=c(9.5,0.5), axes=F, ylab="", xlab="")

axis(1, at=c(1:8), labels=FALSE, las=2)
text(1:8, rep(10.1, 8), colnames(coefs.all_models), xpd=TRUE, srt=90, cex=0.6, col="black")

axis(2, at=c(1:9), cex.lab=0.45, cex.axis=0.45, labels=c(expression("NC"[-1]), expression("MAX"[JUL-1]), expression("MAX"[JUN-1]), expression("MAX"[JUL-2]),
                                                         expression("MAX"[JUN-2]), expression("PRE"[JUL-1]), expression("PRE"[JUN-1]),
                                                         expression("PRE"[JUL-2]), expression("PRE"[JUN-2])), las=1)

mtext("Model coefficents", 2, line=1.9, cex=0.8)
abline(v=(0:15)+0.5, lwd=3.5, col="white")
abline(h=(0:9)+0.5, lwd=3.5, col="white")
box()

mtext(text="B)", side=2,  padj=-8.2, las=1, line=1.5, cex=0.9)

par(new=TRUE)
par(fig=c(0,1,0.0,0.077))

Col<- rbPal(180)
image(x=-25:25, y=1, matrix(1:50), col=Col, ylim=c(8.5,0.5), axes=F, ylab="", xlab="")
axis(1, at=c(-4:4)*10, labels=c(-4:4), tck=0.1, padj=-1)
mtext("Coefficent", 1, line=1, cex=0.7)
box()

dev.off()


### Moving Spearman Rank Correlation ####
# (on longest NC)

# define function for moving correlation with minimum sample size
sample <- 15
my.fun <- function(x,y) {
  my.df <- data.frame(x,y)
  my.df.cmpl <- my.df[complete.cases(my.df), ]
  
  # 15 complete obs is the minimum for cor.test
  if (nrow(my.df.cmpl)<sample) {
    return(rep(NA, 2))
  } else {
    my.test <- cor.test(my.df.cmpl$x,my.df.cmpl$y, method="spearman")
    return(c(my.test$estimate, my.test$p.value))
  }
}

# plot
#tiff("Figure_6_revised.tif", width=18, height=14, units="cm", res=800)
par(mar=c(2,2,0.2,0.2)+0.1, mgp=c(1.1,0.2,0), tck=0.02, lwd=1, cex=1, cex.lab=1, cex.axis=1,  
    pch=16, mfrow=c(3,3))

number <- 28 #window size (28 is the largest window giving 4 independent intervals for the period 1903-2016)
year.start<- 1928 # end of first window
year.end<- 2016 #?

# create masing and climate datasets
mast.MCA<- cbind(mast[,1], mast[,sites])
clim.MCA<- cbind(clim[[2]][,9:10], clim[[2]][,sites]) #temperature
# loop over sites for calculating MCA and creating plot.
for(i in 1:8) {
  x<- mast.MCA[,i+1] # select mast data
  y1<- subset(clim.MCA, clim.MCA[,2] == 19)[,i+2] # select climate data, based on NUTS and the month required
  y2<- subset(clim.MCA, clim.MCA[,2] == 31)[,i+2]
  y3<- subset(clim.MCA, clim.MCA[,2] == 18)[,i+2]
  y4<- subset(clim.MCA, clim.MCA[,2] == 30)[,i+2]
  
  rc.JUL.1<-running(x,y1,fun=my.fun, width=number)#calculate the running correlation
  rc.JUL.1.p<-running(x,y1,fun=my.fun, width=number)#calculate the running correlation
  
  rc.JUL.2<-running(x,y2,fun=my.fun, width=number)#calculate the running correlation
  rc.JUL.2.p<-running(x,y2,fun=my.fun, width=number)#calculate the running correlation
  
  rc.JUN.1<-running(x,y3,fun=my.fun, width=number)#calculate the running correlation
  rc.JUN.1.p<-running(x,y3,fun=my.fun, width=number)#calculate the running correlation
  
  rc.JUN.2<-running(x,y4,fun=my.fun, width=number)#calculate the running correlation
  rc.JUN.2.p<-running(x,y4,fun=my.fun, width=number)#calculate the running correlation
  
  plot(c(1901,2020), c(-0.8,0.8), ylab="Correlation", type="n", xlab="Year", 
       ylim=c(-0.8,0.8), lwd=3, col="red")
  legend("topleft", names(mast.MCA)[i+1], bty="n", xjust=0, inset=0.002)
  abline(h=0)
  
  # plot the result for each month, including plotting the individual segments with line width depednent on significance
  lines(c(1901, 1929, 1957, 1985, 2013), rc.JUN.1[1, c(4, 32, 58, 86, 86)], type="s", lwd=1, col="pink")
  lines(c(1901, 1928), rc.JUN.1[1, c(4, 4)], type= "s", lwd=ifelse(rc.JUN.1[2,4]<0.05,3,1), col="pink", lend=3)
  lines(c(1928, 1957), rc.JUN.1[1, c(32, 32)], type= "s", lwd=ifelse(rc.JUN.1[2,32]<0.05,3,1), col="pink", lend=3)
  lines(c(1957, 1985), rc.JUN.1[1, c(58, 58)], type= "s", lwd=ifelse(rc.JUN.1[2,58]<0.05,3,1), col="pink", lend=3)
  lines(c(1985, 2013), rc.JUN.1[1, c(86, 86)], type= "s", lwd=ifelse(rc.JUN.1[2,86]<0.05,3,1), col="pink", lend=3)
  
  lines(c(1901, 1929, 1957, 1985, 2013), rc.JUN.2[1, c(4, 32, 58, 86, 86)], type="s", lwd=1, col="lightblue")
  lines(c(1901, 1929), rc.JUN.2[1, c(4, 4)], type= "s", lwd=ifelse(rc.JUN.2[2,4]<0.05,3,1), col="lightblue", lend=3)
  lines(c(1929, 1957), rc.JUN.2[1, c(32, 32)], type= "s", lwd=ifelse(rc.JUN.2[2,32]<0.05,3,1), col="lightblue", lend=3)
  lines(c(1957, 1985), rc.JUN.2[1, c(58, 58)], type= "s", lwd=ifelse(rc.JUN.2[2,58]<0.05,3,1), col="lightblue", lend=3)
  lines(c(1985, 2013), rc.JUN.2[1, c(86, 86)], type= "s", lwd=ifelse(rc.JUN.2[2,86]<0.05,3,1), col="lightblue", lend=3)
  
  lines(c(1901, 1929, 1957, 1985, 2013), rc.JUL.2[1, c(4, 32, 58, 86, 86)], type="s", lwd=1, col="blue")
  lines(c(1901, 1929), rc.JUL.2[1, c(4, 4)], type= "s", lwd=ifelse(rc.JUL.2[2,4]<0.05,3,1), col="blue", lend=3)
  lines(c(1929, 1957), rc.JUL.2[1, c(32, 32)], type= "s", lwd=ifelse(rc.JUL.2[2,32]<0.05,3,1), col="blue", lend=3)
  lines(c(1957, 1985), rc.JUL.2[1, c(58, 58)], type= "s", lwd=ifelse(rc.JUL.2[2,58]<0.05,3,1), col="blue", lend=3)
  lines(c(1985, 2013), rc.JUL.2[1, c(86, 86)], type= "s", lwd=ifelse(rc.JUL.2[2,86]<0.05,3,1), col="blue", lend=3)
  
  lines(c(1901, 1929, 1957, 1985, 2013), rc.JUL.1[1, c(4, 32, 58, 86, 86)], type="s", lwd=1, col="red")
  lines(c(1901, 1929), rc.JUL.1[1, c(4, 4)], type= "s", lwd=ifelse(rc.JUL.1[2,4]<0.05,3,1), col="red", lend=3)
  lines(c(1929, 1957), rc.JUL.1[1, c(32, 32)], type= "s", lwd=ifelse(rc.JUL.1[2,32]<0.05,3,1), col="red", lend=3)
  lines(c(1957, 1985), rc.JUL.1[1, c(58, 58)], type= "s", lwd=ifelse(rc.JUL.1[2,58]<0.05,3,1), col="red", lend=3)
  lines(c(1985, 2013), rc.JUL.1[1, c(86, 86)], type= "s", lwd=ifelse(rc.JUL.1[2,86]<0.05,3,1), col="red", lend=3)
  
}


plot(1, 1, ylab="", type="n", xlab="", axes=F)
legend("left", c("MAX.JUL-1", "MAX.JUN-1", "MAX.JUL-2", "MAX.JUN-2", "", "p < 0.05"), col=c("red", "pink", "blue", "lightblue", "White", "black"),
       lwd=c(1,1,1,1,1,3), bty="n")

dev.off()

# create the violin plot
# some plotting parameters
par(mfcol=c(1,1))
b.col1<- "darkred"
b.col2<- "darkblue"
col1<- rgb(255, 15, 15, 50, maxColor=255)
col2<- rgb(15, 15, 255, 50, maxColor=255)

col3<- "pink"
col4<- "lightblue"

r.col1<- col1
r.col2<- col2
m.col1<- b.col1
m.col2<- b.col2
PCH<- 16
xx<- 0.7
xx2<- 1.3
wex= 1
cex.point<- 0.3
fact<- 0.3
point.col1<- "red"
point.col2<- rgb(255, 0, 0, 90, maxColor=255)
point.col3<- "blue"
point.col4<- rgb(0, 0, 255, 90, maxColor=255)

# plot as violin plots
#png("PREC.MCA.28year window.VIOLIN.png", width=8.6, height=6, units="cm", res=800)

par(mar=c(1.5,1.5,0.2,0.2)+0.1, mgp=c(0.6,0.2,0), tck=0.02, lwd=1, cex=0.9, cex.lab=1, cex.axis=1,  
    pch=16)

boxplot((1:4), na.rm=TRUE, border="white", xlim=c(0.5,8.5), ylim=c(-0.81, 0.81), axes=FALSE)
abline(h=0, col="black")
box()
axis(2, cex.axis=0.5, at=c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8), las=1)
axis(1, cex.axis=0.5, padj=-1.2, at=c(1:8), labels=c("DE1", "DE2", "DE9", "DEF", "DK0", "NL1", "SE2", "UKJ"))
mtext("NUTS REGION", 1, line=0.6, cex=0.6)
mtext("Correlation coefficent", 2, line=0.8, cex=0.6)

for(i in 1:8) {
  x<- mast.MCA[,i+1]
  y1<- subset(clim.MCA, clim.MCA[,2] == 19)[,i+2]
  y2<- subset(clim.MCA, clim.MCA[,2] == 31)[,i+2]
  
  rc.JUL.1<-running(x,y1,fun=my.fun, width=number)#calculate the running correlation
  rc.JUL.2<-running(x,y2,fun=my.fun, width=number)#calculate the running correlation  
  
  vioplot(rc.JUL.1[1,][!is.na(rc.JUL.1[1,])], na.rm=TRUE, ylim=c(-0.9,0.9), add=TRUE, at=i, col=col1, border=col3, drawRect=F, wex=wex)
  points(jitter(rep(i, ncol (rc.JUL.1)), amount=fact), rc.JUL.1[1,], cex=cex.point, col=ifelse(rc.JUL.1[2,]<0.05, point.col1, point.col2))
  points(i, median(rc.JUL.1[1,][!is.na(rc.JUL.1[1,])]), pch=PCH, col=b.col1)
  
  vioplot(rc.JUL.2[1,][!is.na(rc.JUL.2[1,])], na.rm=TRUE, ylim=c(-0.9,0.9), add=TRUE, at=i, col=col2, border=col4, drawRect=F, wex=wex)
  points(jitter(rep(i, ncol (rc.JUL.1)), amount=fact), rc.JUL.2[1,], cex=cex.point, col=ifelse(rc.JUL.2[2,]<0.05, point.col3, point.col4))
  points(i, median(rc.JUL.2[1,][!is.na(rc.JUL.2[1,])]), pch=PCH, col=b.col2)
}

###  temporal trends in sensitivity ####
# by testing for interaction between climate variables and Year (in a logistic regression model)

# create data for model

interaction.model<- matrix(nrow=8, ncol=8, byrow=TRUE)
for(i in 1:8) { 
  site<- sites[i]
  data<- data.frame (cbind (mast.ord [,site], MAST.L, TEMP.JUL.1, TEMP.JUN.1, TEMP.JUL.2, TEMP.JUN.2, PREC.JUL.1, PREC.JUN.1, PREC.JUL.2, PREC.JUN.2))
  data$Year<- c(1903:2016)
  
  fit.jun1 <- lrm(data[,1] ~ MAST.L + TEMP.JUN.1*scale(Year), data=data)
  fit.jul1 <- lrm(data[,1] ~ MAST.L + TEMP.JUL.1*scale(Year), data=data)
  fit.jun2 <- lrm(data[,1] ~ MAST.L + TEMP.JUN.2*scale(Year), data=data)
  fit.jul2 <- lrm(data[,1] ~ MAST.L + TEMP.JUL.2*scale(Year), data=data)
  
  interaction.model[i,1]<- data.frame(anova(fit.jun1))[6,3]
  interaction.model[i,2]<- data.frame(anova(fit.jul1))[6,3]
  interaction.model[i,3]<- data.frame(anova(fit.jun2))[6,3]
  interaction.model[i,4]<- data.frame(anova(fit.jul2))[6,3]
  
  interaction.model[i,5]<- fit.jun1$coefficients[8]
  interaction.model[i,6]<- fit.jul1$coefficients[8]
  interaction.model[i,7]<- fit.jun2$coefficients[8]
  interaction.model[i,8]<- fit.jul2$coefficients[8]
}

colnames(interaction.model)<- c("p.jun1", "p.jul1", "p.jun2", "p.jul2", "coef.jun1", "coef.jul1", "coef.jun2", "coef.jul2")
rownames(interaction.model)<- sites

year.model<- matrix(nrow=8, ncol=8, byrow=TRUE)
for(i in 1:8) { 
  site<- sites[i]
  data<- data.frame (cbind (mast.ord [,site], MAST.L, TEMP.JUL.1, TEMP.JUN.1, TEMP.JUL.2, TEMP.JUN.2, PREC.JUL.1, PREC.JUN.1, PREC.JUL.2, PREC.JUN.2))
  data$Year<- c(1903:2016)
  
  fit.jun1 <- lrm(data[,1] ~ MAST.L + TEMP.JUN.1*scale(Year), data=data)
  fit.jul1 <- lrm(data[,1] ~ MAST.L + TEMP.JUL.1*scale(Year), data=data)
  fit.jun2 <- lrm(data[,1] ~ MAST.L + TEMP.JUN.2*scale(Year), data=data)
  fit.jul2 <- lrm(data[,1] ~ MAST.L + TEMP.JUL.2*scale(Year), data=data)
  
  year.model[i,1]<- data.frame(anova(fit.jun1))[4,3]
  year.model[i,2]<- data.frame(anova(fit.jul1))[4,3]
  year.model[i,3]<- data.frame(anova(fit.jun2))[4,3]
  year.model[i,4]<- data.frame(anova(fit.jul2))[4,3]
  
  year.model[i,5]<- fit.jun1$coefficients[7]
  year.model[i,6]<- fit.jul1$coefficients[7]
  year.model[i,7]<- fit.jun2$coefficients[7]
  year.model[i,8]<- fit.jul2$coefficients[7]
}

colnames(year.model)<- c("p.jun1", "p.jul1", "p.jun2", "p.jul2", "coef.jun1", "coef.jul1", "coef.jun2", "coef.jul2")
rownames(year.model)<- sites

### # END OF CODE ####

