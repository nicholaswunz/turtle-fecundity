# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 13/05/2020
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: Turtle reproductive output
# Description: Spaital disribution of nesting sites map

# Load packages into the workspace
install.packages("broom")
library(rgdal) # Bindings for the 'Geospatial' Data Abstraction Library e.g. readOGR(), spTransform()
library(sf) # encoding spatial vector data
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # aligning  graphs
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
turtle.dist <- subset(turtle.dat, !is.na(lon) & !is.na(lat)) # remove studies with no gps points

levels(turtle.dist$red.list)
turtle.dist$red.list <- factor(turtle.dist$red.list, levels = c("Critically endangered", "Endangered", "Vulnerable","Data deficient"))
iucn <- rbind(c("#D6302D","#D1633A","#CA9D00","#6E6E6E"))

# check colour-blind friendly palettes
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
mycol <- brewer.pal(7, "Dark2")
mycol <- c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

# plot with mean mass (continuous) - For Fig S1
ggplot(turtle.dist, aes(x = lat, y = pred.mass, colour = species)) + mytheme +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  xlim(-100, 100) +
  scale_color_manual(values = mycol) + coord_flip() +
  facet_wrap(~ species, ncol = 2, scales = "free_x")

## MPA DATASET ##---------------------------------------------------------------------------------------------------
# upload MPAs using wdpar package
#install.packages("wdpar", repos = "https://cran.rstudio.com/")
library(wdpar)
wdpa.dat <- wdpa_fetch("global", wait = TRUE) # download global protected area data

class(wdpa.dat) # check class type
st_crs(coastlines) # coordinate system
names(wdpa.dat)

# subset marine only data
mpa.raw <- subset(wdpa.dat, MARINE == "2") # subset marine data only
mpa.raw <- subset(mpa.raw, STATUS != "Proposed") # only keep current MPAs
mpa.raw <- subset(mpa.raw, VERIF != "Not Reported") # only include state and expert verified data
mpa.raw <- subset(mpa.raw, ISO3 != "ABNJ") # remove areas beyond national jurisdiction (ABNJ)
names(mpa.raw)

# subset no take zone
notake.raw <- subset(mpa.raw, NO_TAKE %in% c("All","Part"))

# notake.dat <- wdpa_clean(notake.raw) # clean dataset 
# Error: Can't find column `geometry` in `.data`.

attr(notake.raw, "sf_corlumn")

# function to show different classes
types <- vapply(sf::st_geometry(notake.raw), function(x) {
  class(x)[2]}, "")
unique(types) # show class types

# extract MULTIPOLYGON from raw data
polys <- notake.raw[ grepl("*MULTIPOLYGON", types), ]

notake.sp <- as(polys, "Spatial") # convert to sp class
notake.sp2 <- SpatialPolygonsDataFrame(notake.sp, notake.sp@data) #turn the data into a spatial data frame

notake.buf <- st_buffer(polys, dist = 0.05) # 5 km buffer zone
notake.5km <- as(notake.buf, "Spatial") # convert to sp class
#plot(notake.buf)
#plot(notake.sp, add = TRUE,col = "red")

## MAP DATASET ##---------------------------------------------------------------------------------------------------
# download the data
#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip",
#destfile = 'coastlines.zip')

# For coastline only shape file
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
coastlines2 <- SpatialLinesDataFrame(coastlines, coastlines@data) #turn the data into a spatial data frame

world <- readOGR("ne_50m_land/ne_50m_land.shp")
world2 <- SpatialPolygonsDataFrame(world, world@data) #turn the data into a spatial data frame

# check classes of all variables in spatial dataset
sapply(world2@data, class)

# nesting distribution by species (Fig S1)
ggplot() + mytheme + coord_fixed() +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group), colour = "#64686b", 
               fill = "#F1EEE7", size = 0.2) +
  geom_point(data = turtle.dist, aes(x = lon, y = lat, colour = species, fill = species),
             size = 2, alpha = 0.8) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"),
        legend.position = "bottom") +
  facet_wrap(~ species, ncol = 2)

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

table(turtle_in_mpa$species) # how many nesting turtles in MPA's
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
figS2a <- ggplot() + mytheme + coord_fixed() +
  geom_polygon(data = notake.sp2, aes(x = long, y = lat, group = group), colour = "#6baed6", # MPA data
               fill = "#6baed6") +
  geom_polygon(data = world2, aes(x = long, y = lat, group = group), colour = "#64686b", # world space
               fill = "#F1EEE7", size = 0.2) +
  geom_point(data = turtle.dist, aes(x = lon, y = lat), colour = "#81CA9F", # nest distribution
             size = 3, alpha = 0.8) +
  geom_point(data = turtle_in_mpa.df, aes(x = lon, y = lat),colour = "#14505C", # nest in MPAs
             size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"),
        legend.position = "bottom")

# plot mass over species grouped by MPAs
library(dplyr)
count <- turtle.dist %>% 
  count(species, MPA2, sort = TRUE)

figS2b <- ggplot(turtle.dist, aes(x = species, y = pred.mass, colour = MPA2)) + 
  mytheme +
  #geom_boxplot(width = 0.5, position = position_dodge(0.8)) +
  geom_jitter(size = 3.5, position = position_dodge(0.8), alpha = 0.6) +
  scale_colour_manual(values = c("#81CA9F", "#14505C")) +
  ylab("Predicted body mass (kg)") + xlab(NULL) +
  geom_text(data = count, aes(species, Inf, label = n), vjust = 1, size = 3, 
            position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

plot_grid(figS2a, figS2b, ncol = 1, align = "v", axis = "bt", 
          rel_heights = c(1,0.5), labels = c("A", "B"))

# proportion of nesting females in MPAs
library(scales)
# testing percentage
ggplot(count, aes(x = "", y = n, group = MPA2, fill = MPA2))+
  mytheme +
  geom_bar(stat = "identity", position = position_fill()) +
  coord_polar(theta = "y") +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(y = n/2 + c(0, cumsum(n)[-length(n)]), 
                label = percent(n/100)), size = 5) +
  facet_wrap(~ species) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

# For Fig S2b
ggplot(data = count, aes(x = "", y = n, fill = MPA2)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#81CA9F", "#14505C")) +
  facet_wrap( ~ species) + 
  mytheme + theme(axis.text.x = element_blank())


## PROTECTED BEACHES PER COUNTRY ##-----------------------------------------------
# nesting turtle distribution with MPA with protected beaches per country
library(rnaturalearth) # load map
library(rnaturalearthdata) # load data for map

# create dataframe of countries
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# upload mean nesting mass by country
worldturtle <- read.csv("world turtle.csv") # load BPA levels that match the world.csv (include NA)
str(worldturtle)
merged <- merge(world, worldturtle, by.x = "name", by.y = "name")

# plot global protected beaches levels
ggplot(data = merged) + mytheme +
  geom_polygon(data = notake.sp2, aes(x = long, y = lat, group = group), colour = "#6baed6", # MPA data
               fill = "#6baed6") +
  geom_sf(aes(fill = protected.site), colour = NA, size = 0.2, alpha = 0.8) + # % protected beaches
  geom_path(data = coastlines2, aes(x = long, y = lat, group = group), size = 0.1) + # coastline plot
  geom_point(data = turtle.dist, aes(x = lon, y = lat), colour = "black", # add turtle nesting coordinates
             size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"),
        legend.position = "bottom") +
  scale_fill_continuous_sequential(palette = "BluGrn", na.value = "#F1EEE7")

# relationship between % protected beach and body mass
str(worldturtle)
mycol

library(ggrepel) # for adding labels like geom_text_repel
library(brms) # Bayesian Regression Models
library(dplyr)

# GENERAL STAN SPECS
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores avaliable to use

# Caretta caretta
# create regression model via brm
cc.m <- brm(CC_mass ~ CC. + (CC.|temp), data = worldturtle, 
            chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999), max_treedepth = 15)
summary(cc.m)

plot(cc.m, N = 2, ask = FALSE) # Histograms of posterior samples and trace plots of the intercept

# extract posterior samples
cc.post <- brms::posterior_samples(cc.m)
head(cc.post)

# obtain mean intercept and slope
cc.sum <- cc.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_CC.))

# subset 1000 interations
cc.post <- sample_n(cc.post, 1000) 

cc.plot <- ggplot(worldturtle, aes(x = CC., y = CC_mass)) + mytheme +
  geom_abline(data = cc.post, aes(intercept = b_Intercept, slope = b_CC.), # 1000 iterations
              alpha = 0.02, colour = "#D95F02") +
  geom_abline(data = cc.sum, aes(intercept = intercept, slope = slope), # mean fit
              lwd = 1, colour = "#D95F02", alpha = 0.8) +
  geom_point(size = 3, colour = "#D95F02", alpha = 0.8) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  geom_text_repel(aes(label = name), size = 3) +
  ggtitle("Caretta caretta")

# Chelonia mydas
cm.m <- brm(CM_mass ~ CM. + (CM.|temp), data = worldturtle, 
            chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999), max_treedepth = 15)
summary(cm.m)

# extract posterior samples
cm.post <- brms::posterior_samples(cm.m)

# obtain mean intercept and slope
cm.sum <- cm.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_CM.))

# subset 1000 interations
cm.post <- sample_n(cm.post, 1000) 

cm.plot <- ggplot(worldturtle, aes(x = CM., y = CM_mass)) + mytheme +
  geom_abline(data = cm.post, aes(intercept = b_Intercept, slope = b_CM.), # 1000 iterations
              alpha = 0.02, colour = "#1B9E77") +
  geom_abline(data = cm.sum, aes(intercept = intercept, slope = slope), # mean fit
              lwd = 1, colour = "#1B9E77", alpha = 0.8) +
  geom_point(size = 3, colour = "#1B9E77", alpha = 0.8) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  geom_text_repel(aes(label = name), size = 3) +
  ggtitle("Chelonia mydas")

# Dermochelys coriacea
dc.m <- brm(DC_mass ~ DC. + (DC.|temp), data = worldturtle, 
            chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999), max_treedepth = 15)
summary(dc.m)

# extract posterior samples
dc.post <- brms::posterior_samples(dc.m)
head(dc.post)

# obtain mean intercept and slope
dc.sum <- dc.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_DC.))

# subset 1000 interations
dc.post <- sample_n(dc.post, 1000) 

dc.plot <- ggplot(worldturtle, aes(x = DC., y = DC_mass)) + mytheme +
  geom_abline(data = dc.post, aes(intercept = b_Intercept, slope = b_DC.), # 1000 iterations
              alpha = 0.02, colour = "#7570B3") +
  geom_abline(data = dc.sum, aes(intercept = intercept, slope = slope), # mean fit
              lwd = 1, colour = "#7570B3", alpha = 0.8) +
  geom_point(size = 3, colour = "#7570B3", alpha = 0.8) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  geom_text_repel(aes(label = name), size = 3) +
  ggtitle("Dermochelys coriacea")

# Eretmochelys imbricata
# create regression model via brm
ei.m <- brm(EI_mass ~ EI. + (EI.|temp), data = worldturtle, 
            chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999), max_treedepth = 15)
summary(ei.m)

# extract posterior samples
ei.post <- brms::posterior_samples(ei.m)
head(ei.post)

# obtain mean intercept and slope
ei.sum <- ei.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_EI.))

# subset 1000 interations
ei.post <- sample_n(ei.post, 1000) 

ei.plot <- ggplot(worldturtle, aes(x = EI., y = EI_mass)) + mytheme +
  geom_abline(data = ei.post, aes(intercept = b_Intercept, slope = b_EI.), # 1000 iterations
              alpha = 0.02, colour = "#E7298A") +
  geom_abline(data = ei.sum, aes(intercept = intercept, slope = slope), # mean fit
              size = 1, colour = "#E7298A", alpha = 0.8) +
  geom_point(size = 3, colour = "#E7298A", alpha = 0.8) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  geom_text_repel(aes(label = name), size = 3) +
  ggtitle("Eretmochelys imbricata")

# Lepidochelys kempii
lk.plot <- ggplot(worldturtle, aes(x = LK., y = LK_mass)) + mytheme +
  geom_point(size = 3, colour = "#66A61E") +
  #geom_smooth(method = "lm", colour = "#66A61E", size = 1) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  geom_text_repel(aes(label = name), size = 3) +
  ggtitle("Lepidochelys kempii")

# Lepidochelys olivacea
# create regression model via brm
lo.m <- brm(LO_mass ~ LO. + (LO.|temp), data = worldturtle, 
            chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999), max_treedepth = 15)
summary(lo.m)

library(bayesplot)
library(tidyverse)
worldturtle2 <- worldturtle %>% drop_na(LO_mass, LO.)
y      <- worldturtle2$LO_mass
yrep   <- posterior_predict(lo.m, draws = 500)
ppc_scatter_avg(y, yrep, alpha = 0.7)

# extract posterior samples
lo.post <- brms::posterior_samples(lo.m)
head(lo.post)

# obtain mean intercept and slope
lo.sum <- lo.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_LO.))

# subset 1000 interations
lo.post <- sample_n(lo.post, 1000) 

lo.plot <- ggplot(worldturtle, aes(x = LO., y = LO_mass)) + mytheme +
  geom_abline(data = lo.post, aes(intercept = b_Intercept, slope = b_LO.), # 1000 iterations
              alpha = 0.02, colour = "#E6AB02") +
  geom_abline(data = lo.sum, aes(intercept = intercept, slope = slope), # mean fit
              size = 1, colour = "#E6AB02", alpha = 0.8) +
  geom_point(size = 3, colour = "#E6AB02", alpha = 0.8) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  geom_text_repel(aes(label = name), size = 3) +
  ggtitle("Lepidochelys olivacea")

# Natator depressa (only 1 datapoint)

# combine all species plot
plot_grid(cc.plot, cm.plot, dc.plot, ei.plot, lk.plot, lo.plot,
          ncol = 2, align = "hv", axis = "tblr")

