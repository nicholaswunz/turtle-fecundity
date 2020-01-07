# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 28/11/2019
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: Turtle reproductive output
# Description: Data analysis & figure production

# install and load packages
install.packages("")
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # organise figures
library(plyr) # for ddply summary stats
library(colorspace)
# detach("package:dplyr", unload=TRUE) # clashes with dplyr at times
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

## ACROSS SPECIES ANALYSIS ##-------------------------------------------------
# load data
turtle.dat <- read.csv("reproductive output.csv")
str(turtle.dat)
names(turtle.dat)
length(unique(turtle.dat$egg.volume)) # check number of data with egg volume

# following code based on Barneche et al 2018 Si eq 2 and 3

# total clutch volume = fecundity * egg-volume
# total clutch volume = b0_m1 * mass ^ b1_m1 * b0_m2 * mass ^ b1_m2 # written as model
# total clutch volume = b0_m1 * b0_m2 * mass ^ (b1_m1 + b1_m2) # rearranged for better order

b0_m1 * mass ^ b1_m1 # fecundity
b0_m2 * mass ^ b1_m2 # egg-volume

# obtain predicted mean mass from each species
spp.mass <- aggregate(turtle.dat$pred.mass, list(turtle.dat$species), mean, na.rm = TRUE)

## For overall clutch volume ##
# calculate fercundity slope - all species *pred.clutch = predicted total clutch size
levels(turtle.dat$species) # check species list
fecundity.m  <-  lm(log(pred.clutch) ~ log(pred.mass), data = turtle.dat)
summary(fecundity.m)

exp(coef(fecundity.m)['(Intercept)']) # unlog intercept
coef(fecundity.m)['log(pred.mass)'] # unlog slope

# caculate egg volume slope - all species
eggvolume.m  <-  lm(log(egg.volume) ~ log(pred.mass), data = turtle.dat)
summary(eggvolume.m)

# calculate total clutch volume slope from fecundity and egg volume slope
clutchvolume.m <- coef(fecundity.m)['log(pred.mass)'] + 
  coef(eggvolume.m)['log(pred.mass)'] # slope calculation
coef(fecundity.m)['(Intercept)'] * coef(eggvolume.m)['(Intercept)'] #intercept calculation

## Caluclate clutch volume for each species ##
# Calculate fecundity slope for each species
# Break up turtle.dat by species, then fit the specified model to each piece & return a list
fecundity.m <- dlply(turtle.dat, "species", function(df)
  lm(log(pred.clutch) ~ log(pred.mass), data = df))

# Apply coef to each model and return a data frame
fecundity.m2 <- ldply(fecundity.m, coef)
fecundity.m2 <- rename(fecundity.m2, c("log(pred.mass)" = "fecundity.slope"))

# Calculate egg volume slope for each species
eggvolume.m <- dlply(turtle.dat, "species", function(df)
  lm(log(egg.volume) ~ log(pred.mass), data = df))

eggvolume.m2 <- ldply(eggvolume.m, coef)
eggvolume.m2 <- rename(eggvolume.m2, c("log(pred.mass)" = "eggvolume.slope"))

slope.m <- merge(fecundity.m2, eggvolume.m2, by = "species") # merge both model outcomes
slope.m$clutchvolume.int <- slope.m$"(Intercept).x" * slope.m$"(Intercept).y"
slope.m$clutchvolume.slope <- slope.m$fecundity.slope +
  slope.m$eggvolume.slope # calculate clutch volume slope for each specoes


library(dplyr)
# calculate mean predicted mass for each species
mean.mass <- turtle.dat %>%
  group_by(species) %>%
  summarize(mean = mean(pred.mass, na.rm = TRUE))

data.frame(mean.mass)

# subset mean length to slope.m
slope.m$mass <- mean.mass$mean
row.names(slope.m) <- slope.m$species

# if total energy clutch slope is >1, then indicates reproductive hyperallometry

str(slope.m)
# plot clutch volume with species (categorical)
ggplot(slope.m, aes(x = species, y = clutchvolume.slope, colour = clutchvolume.slope)) + mytheme +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3)

# plot clutch volume with mean predicted mass (continuous)
ggplot(slope.m, aes(x = log(mass), y = clutchvolume.slope, colour = species)) + mytheme +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3) + geom_smooth(method='lm')

# plot female length vs clutch

mass.plot <- ggplot(turtle.dat, aes(x = log10(pred.mass), y = log10(pred.clutch), fill = species, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.3, show.legend = FALSE) +
  geom_smooth(aes(fill = species), method = lm, se = FALSE, size = 2) +
  xlab(expression("Predicted body mass (10"^"x"*", kg)")) + 
  ylab(expression("Predicted fecundity (total clutch output, 10"^"y"*")")) +
  scale_fill_discrete_sequential(palette = "BluGrn", nmax = 7) +
  scale_colour_discrete_sequential(palette = "BluGrn", nmax = 7)

  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)
  
# plot female mass vs eggs size
volume.plot <- ggplot(turtle.dat, aes(x = log10(pred.mass), y = log10(egg.volume), fill = species, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.3, show.legend = FALSE) +
  geom_smooth(aes(fill = species), method = lm, se = FALSE, size = 2) +
  xlab(expression("Predicted body mass (10"^"x"*", kg)")) + 
  ylab(expression("Egg volume (10"^"y"*", cm"^"3"*")")) +
  scale_fill_discrete_sequential(palette = "BluGrn", nmax = 7) +
  scale_colour_discrete_sequential(palette = "BluGrn", nmax = 7)

prow <- plot_grid(
  mass.plot + theme(legend.position="none"),
  volume.plot + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)

# create legend from mass.plot to plot under the graph
legend_b <- get_legend(mass.plot +
                         guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))

# add the legend to the row made earlier. Give it one-tenth of the height of one plot (via rel_heights).
output.plot <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

## PHYLOGENETIC RECONSTRUCTION ##-------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")
library(rotl) # Interface to the 'Open Tree of Life' API.
library(ape) # deal with phylogenetic data
library(ggtree) # plotting trees in ggplot2 format
library(phytools) # merge data with tree and phylogenetic signal calculation
library(ggstance) # associate graphs to phylogenetic trees

slope.m$species <- as.character(slope.m$species) # sort taxa names as character first

# generating list of species
species <- sort(unique(slope.m$species))

# obtaining dataframe listing the Open Tree identifiers potentially matching our list of species.
taxa <- tnrs_match_names(names = species)
taxa

# check if species list matcg OT identifier
taxa[taxa$approximate_match == TRUE,] # none

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = ott_id(taxa),label_format = "name")

plot(tree, cex = 0.6, label.offset = 0.1, no.margin = TRUE) # plot tree

# If polytomies exist, the output will be `FALSE`, and vice versa.
is.binary.tree(tree)

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# convert tree to correlation matrix
phylo_cor <- vcv(phylo_branch, cor = T)
plot(phylo_branch)
tree <- phylo_branch

save(phylo_cor, file = "C:/Users/chwu6540/Dropbox (Personal)/Turtle project/Data/phylo_cor_all.Rdata")


setdiff(slope.m$species, as.character(tree$tip.label)) # listed in the database but not in the tree
setdiff(as.character(tree$tip.label),slope.m$species) # listed in the tree but not in the database
# all ok!

both <- intersect(tree$tip.label, rownames(slope.m)) # matching tree tip label with rownames
tree2 <- drop.tip(tree, setdiff(both, tree$tip.label)) # drop any tips not matching rownames
data <- slope.m[both,] # attach intersect with dataframe

slope <- as.matrix(slope.m)[,7] # extract trait (clutch volume slope)
mass <- as.matrix(slope.m)[,8] # extract trait (predicted mass)
fit <- phytools::fastAnc(tree2, slope, vars = TRUE, CI = TRUE) # calculate fast estimation of ML ancestral states

td <- data.frame(node = nodeid(tree2, names(slope)),
                 trait = slope) # combine length with tree nodes
td$trait <- as.numeric(levels(td$trait))[td$trait] # change to numeric
nd <- data.frame(node = names(fit$ace), trait = fit$ace) # combine ancestral states with tree nodes

d <- rbind(td, nd)
d$node <- as.numeric(d$node) # change to numeric
tree3 <- full_join(tree2, d, by = 'node') # join trait with tree

# plot tree
tree.plot <- ggtree(tree3, aes(color = trait), ladderize = FALSE, size = 2) +
  geom_tiplab(size = 3, hjust = -0.1) +
  scale_colour_continuous_sequential(palette = "BluGrn")

sequential_hcl(7, "BluGrn") #testing

# combine with slope points
slope.plot <- facet_plot(tree.plot + xlim_tree(1.3), panel = "Clutch volume slope",
           data = data, geom = geom_point, mapping = aes(x = clutchvolume.slope), size = 4) +
  theme_tree2() + mytheme +
  geom_vline(xintercept = 1, linetype = "dashed")

# phylogenetic signal
trait <- slope.m[,7] # extract clutch volume slope
names(trait) <- slope.m[,1] # extract species name

# test Pagels lambda with 999 randomizations:
phylosig(phylo_branch, trait, method = "lambda", test = TRUE, nsim = 999)

plot_grid(output.plot, slope.plot, labels = c("", "C"), nrow = 2)


## WITHIN C. MYDAS POPULATION ##----------------------------------------------------
# load dataset
nest.dat <- read.csv("nesting year.csv")
str(nest.dat)
names(nest.dat)

library(ggpointdensity) # plot density in ggplot2
library(quantreg) # plot quartile regression

nest.plot <- ggplot(nest.dat, aes(x = mass, y = sum_n.of.eggs)) + mytheme +
  geom_pointdensity(adjust = 4, size = 3) +
  stat_smooth(method = "lm", size = 2, col = "#14505C", fill = "#C7E5BE") +
  #geom_quantile(quantiles = c(0.5, 0.95), size = 2, col = "#14505C", fill = "#C7E5BE") +
  ylab(expression(paste("Fecundity (total clutch yr"^"-1"*")"))) + 
  xlab("Female mass, Mi, (kg)") +
  scale_colour_continuous_sequential(palette = "BluGrn")


# calculate fercundity slope
total.nest.m <- lm(log10(sum_n.of.eggs) ~ log10(mass), data = nest.dat)
summary(total.nest.m)
total.nest.m$coefficients

10^coef(total.nest.m)['(Intercept)'] # convert to power function

# calculate fercundity slope 95%
#total.nest.m95 <- rq(log10(nest.dat$mean_n.of.eggs) ~ log10(nest.dat$mass), tau = c(0.50,0.95))
#summary(total.nest.m95)
#total.nest.m95$coefficients

# calculate egg volume slope - from species dataset
eggvolume.m <- lm(log10(egg.volume) ~ log10(mass), 
                    data = turtle.dat, species == "Chelonia_mydas")
summary(eggvolume.m)

# calculate egg energy slope
energy.dat <- read.csv("egg energy.csv") # open egg energy file
str(energy.dat)
eggenergy.m <- lm(log10(egg.energy) ~ log10(egg.mass), 
                  data = energy.dat, species == "Chelonia_mydas")
summary(eggenergy.m)

# predict expected egg energy for a female of given mass (check with Diego if this is correct)
# create new empty data.frame
eggenergy.m2 <- data.frame(matrix(ncol = 2, nrow = 1))
colnames(eggenergy.m2) <- c("intercept", "slope") # create column names
eggenergy.m2$intercept <- coef(eggvolume.m)['(Intercept)'] * 
  coef(eggenergy.m)['(Intercept)'] # calculate intercept
eggenergy.m2$slope <- coef(eggvolume.m)['log10(mass)'] * 
  coef(eggenergy.m)['log10(egg.mass)'] # calculate slope

# calculate total energy output
greenint.m <- coef(total.nest.m)['(Intercept)'] * eggenergy.m2$intercept
greenint.m2 <- 10^greenint.m # convert to power function
greenslope.m <- coef(total.nest.m)['log10(mass)'] + eggenergy.m2$slope

# calculate predicted isometry and allometry relationship for C. mydas
mass <-  seq(min(nest.dat$mass - 10), max(nest.dat$mass + 10)) # in cm
isometry <- (greenint.m2 * mass ^ 1) # predicted if isometric, also converted kJ
hyperallometry <- (greenint.m2 * mass ^ greenslope.m) # from results

# plot predicted reproductive-energy output (similar to Diego's Fig 1)
pred.plot <- ggplot() + mytheme +
  geom_smooth(aes(x = mass, y = isometry), size = 1.5, colour = "black", linetype = "dashed") +
  geom_smooth(aes(x = mass, y = hyperallometry), size = 2, colour = "#348781") +
  xlab("Female mass, Mi, (kg)") + ylab("Reproductive output, Ri, (MJ)")

plot_grid(nest.plot, pred.plot, align = "v", axis = "lr", labels = c("A", "B"), ncol = 2)


## TEMPORAL VARIATION ##--------------------------------------------------------------------------
# changes in mass mass per year
library(plyr) # obtain summary stats
nest.mass <- ddply(nest.dat, "year", summarise,
                N = length(mass),
                mean = mean(mass),
                sd = sd(mass),
                se = sd/sqrt(N))
names(nest.mass)

yearmass.plot <- ggplot(nest.mass, aes(x = year, y = mean)) + mytheme +
  geom_smooth(method = lm, size = 2, col = "#14505C", fill = "#C7E5BE", alpha = 0.4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, colour = year), width = 0.2, size = 1) +
  geom_line(linetype = "dashed") +
  geom_point(aes(colour = year), size = 3.5) + xlim(1992, 2020) +
  ylab("Female mass (kg)") + xlab("Year") +
  scale_colour_continuous_sequential(palette = "BluGrn")

str(nest.dat)
yearmass.m <- lm(mass ~ year, data = nest.dat)
summary(yearmass.m)
predict(yearmass.m)

# n of eggs vs year
nest.egg <- ddply(nest.dat, "year", summarise,
             N = length(mean_n.of.eggs),
             mean = mean(mean_n.of.eggs),
             sd = sd(mean_n.of.eggs),
             se = sd/sqrt(N))

yearegg.plot <- ggplot(nest.egg, aes(x = year, y = mean)) + mytheme +
  geom_smooth(method = lm, size = 2, col = "#14505C", fill = "#C7E5BE", alpha = 0.4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, colour = year), width = 0.2, size = 1) +
  geom_line(linetype = "dashed") +
  geom_point(aes(colour = year), size = 3.5) + xlim(1992, 2020) +
  ylab("Mean clutch size") + xlab("Year") +
  scale_colour_continuous_sequential(palette = "BluGrn")

yearegg.m <- lm(mean_n.of.eggs ~ year, data = nest.dat)
summary(yearegg.m)

# group plots
plot_grid(yearmass.plot, yearegg.plot, labels = c('A','B'), nrow = 2, align = "h")


# plot changes in body length per year (distribution plot)
install.packages("ggridges")
library(ggridges)

nest.dat2 <- transform(nest.dat, year2 = as.factor(year)) # for geom_density to work

ggplot(nest.dat2, aes(x = mass, y = year2, fill = year2)) + mytheme +
  geom_density_ridges(scale = 1, rel_min_height = 0.01) +
  xlab("Female mass (kg)") + ylab("Year") +
  scale_fill_discrete_sequential(palette = "BluGrn")



## EXTRA/TEST CODES ##-------------------------------------------------------

# example code from Diego
b0_m1  <-  0.5 # incercept
b1_m1  <-  0.9 # mass or length slope
length <-  seq(80, 120) # in cm
mass   <-  5 * length ^ 3 # estimated mass based on square-cube law (mass = length ^3)
plot(mass ~ length) # plot relationship to show non-linear relationship

fecundity  <-  b0_m1 * mass ^ b1_m1 # fecundity
log(b0_m1) + b1_m1 * log(mass) # log slope and mass

# run model
model_fecundity  <-  lm(log(fecundity) ~ log(mass))
summary(model_fecundity)

exp(coef(model_fecundity)['(Intercept)']) # unlog intercept
coef(model_fecundity)['log(mass)'] # unlog slope

# different type of plot
tree.plot2 <- ggtree(tree2, size = 2, layout = "circular")
data1 <- data[c(1)] # categorical traits
data2 <- data[c(-1,-2,-4,-7)] # continuous traits
# plot for categorical
p1 <- gheatmap(tree.plot2, data1, width = 0.2,
               colnames_angle = 95, colnames_offset_y = 0.25) +
  scale_fill_brewer("Reds")

library(ggnewscale)
# plot for continuous
p2 <- p1 + new_scale_fill()
gheatmap(p2, data2, offset = 0.3, width = 0.3,
         colnames_angle = 90, colnames_offset_y = 0.25) +
  scale_fill_distiller("Blues", direction = 1)
