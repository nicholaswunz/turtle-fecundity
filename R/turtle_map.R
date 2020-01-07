# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 16/11/2019
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: Turtle reproductive output
# Description: Dreating nesting distribution map

# Load packages into the workspace
install.packages("broom")
library(rgdal) # Bindings for the 'Geospatial' Data Abstraction Library e.g. readOGR(), spTransform()
library(sf) # encoding spatial vector data
library(ggplot2) # use ggplot2 method of plotting
library(cowplot)
mytheme <- theme_bw() + {theme(panel.border = element_blank(), # Border around plotting area.
                               panel.grid.major = element_blank(), # Major grid lines blank
                               panel.grid.minor = element_blank(), # Minor grid lines blank
                               axis.line = element_line(colour = "black", size = 0.8), # axis lin size
                               axis.ticks = element_line(colour = "black", size = 0.8),
                               axis.text = element_text(size = 10, colour = "black"), # axis text size
                               axis.title=element_text(size = 10), #axis title size
                               panel.background = element_rect(fill = "transparent"), # bg of the panel
                               plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                               legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                               legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg)
} # set up plot theme

# set directory
setwd('C:/Users/nicho/Dropbox (Personal)/Turtle project/Data') # on laptop
setwd('C:/Users/chwu6540/Dropbox (Personal)/Turtle project/Data') # on desktop

## TURTLE DATASET ##------------------------------------------------------------------------------------------------
# load data
turtle.dat <- read.csv("reproductive output.csv")
str(turtle.dat)
names(turtle.dat)
turtle.dist <- subset(turtle.dat, !is.na(lon) & !is.na(lat))

levels(turtle.dist$red.list)
turtle.dist$red.list <- factor(turtle.dist$red.list, levels = c("Critically endangered", "Endangered", "Vulnerable","Data deficient"))
iucn <- rbind(c("#D6302D","#D1633A","#CA9D00","#6E6E6E"))

# plot with mean length (continuous)
ggplot(turtle.dist, aes(x = lat, y = pred.mass, colour = red.list)) + mytheme +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  xlim(-100, 100) +
  scale_color_manual(values = iucn) + coord_flip() +
  facet_wrap(~ species, ncol = 2, scales = "free_x")

## MPA DATASET ##---------------------------------------------------------------------------------------------------
# upload MPAs using wdpar package
#install.packages("wdpar", repos = "https://cran.rstudio.com/")
library(wdpar)
wdpa.dat <- wdpa_fetch("global", wait = TRUE) # download global protected area data

class(notake.raw) # check class type
st_crs(coastlines) # coordinate system
names(wdpa.dat)
list(mpa.raw$ISO3)

# subset marine only data
mpa.raw <- subset(wdpa.dat, MARINE == "2") # subset marine dataonly
mpa.raw <- subset(mpa.raw, STATUS != "Proposed") # only kep current MPAs
mpa.raw <- subset(mpa.raw, VERIF != "Not Reported") # only include state and expert verified data
mpa.raw <- subset(mpa.raw, ISO3 != "ABNJ") # remove areas beyond national jurisdiction (ABNJ)
names(mpa.raw)

# subset no take zone
notake.raw <- subset(mpa.raw, NO_TAKE %in% c("All","Part"))
list(notake.raw$ISO3)

notake.dat <- wdpa_clean(notake.raw) # clean dataset 
# Error: Can't find column `geometry` in `.data`.

attr(notake.raw, "sf_corlumn")

# check classes of all variables in spatial dataset
sapply(world2@data, class)

# function to show different classes
types <- vapply(sf::st_geometry(notake.raw), function(x) {
  class(x)[2]
}, "")
unique(types) # show class types

# extract MULTIPOLYGON from raw data
polys <- notake.raw[ grepl("*MULTIPOLYGON", types), ]

# convert to sp class
notake.sp <- as(polys, "Spatial") 
notake.sp2 <- SpatialPolygonsDataFrame(notake.sp, notake.sp@data) #turn the data into a spatial data frame

notake.buf <- st_buffer(polys, dist = 0.05) # 5 km buffer zone
notake.5km <- as(notake.buf, "Spatial")
plot(notake.buf2)
plot(notake.sp,add=TRUE,col="red")

## MAP DATASET ##---------------------------------------------------------------------------------------------------
# download the data
#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip",
#destfile = 'coastlines.zip')

# For coastline only shape file
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
coastlines2 <- SpatialLinesDataFrame(coastlines, coastlines@data) #turn the data into a spatial data frame

world <- readOGR("ne_50m_land/ne_50m_land.shp")
world2 <- SpatialPolygonsDataFrame(world, world@data) #turn the data into a spatial data frame

# nesting distribution by species
ggplot() + mytheme + coord_fixed() +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group), colour = "#64686b", 
               fill = "#F1EEE7", size = 0.2) +
  geom_point(data = turtle.dist, aes(x = lon, y = lat, colour = red.list),
             size = 2, alpha = 0.8) +
  scale_color_manual(values = iucn, na.value = "#F1EEE7") +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"),
        legend.position = "bottom") +
  facet_wrap(~ species, ncol = 2)


library(colorspace)
# nesting disribution with MPA
ggplot() + mytheme + coord_fixed() +
  geom_polygon(data = notake.sp2, aes(x = long, y = lat, group = group), colour = "#6baed6", 
               fill = "#6baed6") +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group), colour = "#64686b", 
               fill = "#F1EEE7", size = 0.2) +
  geom_point(data = turtle.dist, aes(x = lon, y = lat, colour = species),
             size = 3, alpha = 0.8) +
  scale_colour_discrete_sequential(palette = "BluGrn", nmax = 7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"),
        legend.position = "bottom")

## OVERLAP TURTLE DATA WITH MPA ##------------------------------------------------
# convert dataframe to spatial points matching coordinates
st_crs(notake.raw) # check coordinate system
head(turtle.dist)
str(turtle.dist)

# change turtle data to spatial points
xy <- turtle.dist[,c(11,10)] # in long lat order
turtle.sp <- SpatialPointsDataFrame(coords = xy, data = turtle.dist,
                              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# overlap nesting turtles in MPAs (+ 5km)
turtle.sf <- spTransform(turtle.sp, CRS(proj4string(notake.5km)))
turtle_in_mpa <- turtle.sp[notake.5km, ]
plot(turtle_in_mpa)

# number of turtles in MPA and number species in MPA
ISO.count <- over(turtle.sp, notake.sp)
turtle.count <- over(notake.sp, turtle.sp)
table(ISO.count$ISO3) # how many nesting turtles in MPA's 

table(turtle.count$species) # how many nesting turtles in MPA's
table(turtle.dist$species)


# convert turtle in MPA data back to data.frame
turtle_in_mpa.df <- as.data.frame(turtle_in_mpa)
str(turtle_in_mpa.df)
nrow(turtle_in_mpa.df)
# add to csv sheet and reload reprodcutiveoutput.csv

# turtles within MPA
sequential_hcl(5, palette = "BluGrn")
levels(turtle.dist$MPA)

# plot turtle MPA disribution map
fig1a <- ggplot() + mytheme + coord_fixed() +
  geom_polygon(data = notake.sp2, aes(x = long, y = lat, group = group), colour = "#6baed6", 
               fill = "#6baed6") +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group), colour = "#64686b", 
               fill = "#F1EEE7", size = 0.2) +
  geom_point(data = turtle.dist, aes(x = lon, y = lat), colour = "#81CA9F",
             size = 3, alpha = 0.8) +
  geom_point(data = turtle_in_mpa.df, aes(x = lon, y = lat),colour = "#14505C",
             size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"),
        legend.position = "bottom")

# plot mass over species grouped by MPAs
library(dplyr)
count <- turtle.dist %>% 
  count(species, MPA2, sort = TRUE)

fig1b <- ggplot(turtle.dist, aes(x = species, y = pred.mass, colour = MPA2)) + 
  mytheme +
  geom_boxplot(width = 0.5, position = position_dodge(0.8)) +
  geom_jitter(size = 3.5, position = position_dodge(0.8), alpha = 0.6) +
  scale_colour_manual(values = c("#81CA9F", "#14505C")) +
  ylab("Predicted body mass (kg)") + xlab(NULL) +
  geom_text(data = count, aes(species, Inf, label = n), vjust = 1, size = 3, 
            position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

plot_grid(fig1a, fig1b, ncol = 1, align = "v", axis = "bt", 
          rel_heights = c(1,0.5), labels = c("A", "B"))

# proportion of nesting females in MPAs
library(scales)
ggplot(count, aes(x = "", y = n, group = MPA2, fill = MPA2))+
  mytheme +
  geom_bar(stat = "identity", position = position_fill()) +
  coord_polar(theta = "y") +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(y = n/2 + c(0, cumsum(n)[-length(n)]), 
                label = percent(n/100)), size = 5) +
  facet_wrap(~ species) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

ggplot(data = count, aes(x = "", y = n, fill = MPA2)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#81CA9F", "#14505C")) +
  facet_wrap( ~ species) + 
  mytheme + theme(axis.text.x = element_blank())

## ISO TEST ##--------------------------------------------------------------------------------------------------
ggplot(data = notake.raw, aes(x = ISO3, y = log10(REP_M_AREA))) + 
  geom_jitter(position = position_jitter(0)) +
  mytheme + coord_flip() + 
  stat_summary(fun.data = mean_sdl, mult=1, 
               geom = "pointrange", color = "red")

# sum area and n of MPA per ISO
notake.sum <- notake.raw %>%
  as.data.frame() %>%
  select(-Shape) %>%
  group_by(ISO3) %>%
  summarize(area_km = sum(REP_M_AREA), count = n_distinct(REP_M_AREA)) %>%
  ungroup() %>%
  arrange(desc(area_km))

# plot sum area across ISO (categorical)
ggplot(data = notake.sum, aes(x = reorder(ISO3, area_km), y = log10(area_km))) + 
  mytheme + coord_flip() + 
  geom_jitter(position = position_jitter(0))

# plot MPAs on map - colour by ISO
ggplot() + mytheme +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group)) +
  geom_sf(data = notake.raw, aes(fill = ISO3), alpha = 0.6) + 
  geom_point(data = turtle.dist, aes(x = lon, y = lat),
             size = 3, alpha = 0.8, colour = "red") +
  theme(legend.position = "bottom")

year <- subset(notake.raw, STATUS_YR != 0)

# convert to year
year$STATUS_YR <- as.Date(as.character(year$STATUS_YR), format = "%Y")
ggplot(data = year, aes(x = STATUS_YR)) + mytheme +
  geom_histogram()


