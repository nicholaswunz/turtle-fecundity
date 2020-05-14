# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 08/05/2020
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: Turtle reproductive output
# Description: Data analysis & figure production

# install and load packages
#install.packages("brms")
library(brms)
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # organise figures
#library(plyr) # for ddply summary stats
library(dplyr) # for data manipulation like pipes %>%
library(colorspace)
# detach("package:plyr", unload=TRUE) # clashes with dplyr at times
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
setwd('C:/Users/chwu6540/Dropbox (Personal)/Turtle project/Data') # on Desktop

# load data
turtle.dat <- read.csv("reproductive output.csv")
str(turtle.dat)

# Summary of data
nrow(turtle.dat) # number of effect sizes
length(unique(turtle.dat$paper.ID)) # number of studies
length(unique(turtle.dat$country)) # number of countries

min(turtle.dat$egg.volume[!is.na(turtle.dat$egg.volume)])
max(turtle.dat$egg.volume[!is.na(turtle.dat$egg.volume)])

# parameter counts for species
turtle.dat %>% 
  group_by(species) %>% 
  summarise(CCL = length(length[!is.na(length)]),
            mass = length(mass[!is.na(mass)]),
            egg.diam = length(egg.diameter[!is.na(egg.diameter)]),
            egg.mass = length(egg.mass[!is.na(egg.mass)]),
            hatch.length = length(hatchling.length[!is.na(hatchling.length)]),
            hatch.mass = length(hatchling.mass[!is.na(hatchling.mass)]),
            clutch.size = length(clutch.size[!is.na(clutch.size)]),
            clutch.mass = length(clutch.mass[!is.na(clutch.mass)]),
            clutch.freq = length(clutch.freq[!is.na(clutch.freq)])
            )

## PHYLOGENETIC RECONSTRUCTION ##-------------------------------------------------
library(rotl) # Interface to the 'Open Tree of Life' API.
library(ape) # deal with phylogenetic data
library(ggtree) # plotting trees in ggplot2 format; install.packages("BiocManager"); BiocManager::install("ggtree")
library(phytools) # merge data with tree and phylogenetic signal calculation
library(ggstance) # associate graphs to phylogenetic trees

# sort taxa names as character first
turtle.dat$species <- as.character(turtle.dat$species) 

# generating list of species
species <- sort(unique(turtle.dat$species))

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

## ACROSS SPECIES ANALYSIS ##------------------------------------------------------
# Calculate fecundity slope
# GENERAL STAN SPECS
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores avaliable to use

turtle.dat$lnMass    <-  log(turtle.dat$pred.mass)
turtle.dat$lnClutch  <-  log(turtle.dat$pred.clutch)
turtle.dat$animal    <-  as.factor(turtle.dat$species) # rename animal for phylogenetic random effects, vs within species effect

# Test for phylogenetic effects on the intercept, but also a repeated measures effect (i.e. within-species variance) on both the intercept and slope
fecundity.m <- brm(
  lnClutch ~ lnMass + (1 | animal) + (1 + lnMass | species), data = turtle.dat, 
  family = gaussian(), cov_ranef = list(animal = phylo_cor),
  prior = c(
    prior(normal(0, 10), "b"), # mean of 0 and SD of 10 (wide distribution)
    prior(normal(0, 10), "Intercept"), # mean of 0 and SD of 10 (if log then wide negative disribution becomes infinite)
    prior(student_t(3, 0, 20), "sd"), # class of random effect deviation to calcuate - has to be postiive (no negative SD)
    prior(student_t(3, 0, 20), "sigma") # residual SD parameter
  ),
  chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15)
)
# (1 | animal) Phylogentic covariance matrix intercept - may have effect on intercept 
# (1 + lnMass | species) estimate slope variance within each species 

# chains = number of Markov chains
# iter = number of total iterations per chain
# adapt_delta = decrease the number of divergent transitions if warning message occurs
# max_treedepth >10 when depth of tree evaluated in each iteration is exceeded

summary(fecundity.m) # show summary of results

plot(fecundity.m, N = 2, ask = FALSE) # Histograms of posterior samples and trace plots of the intercept
plot(conditional_effects(fecundity.m), points = TRUE) 

# estimate phylogenetic signal (delta) via hypothesis method
hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(fecundity.m, hyp, class = NULL))
# delta = 0.7 under estimates

get_variables(fecundity.m) # list of obtainable variables

fixefs  <-  fixef(fecundity.m) # fixed effects summary (all data)
ranefs  <-  ranef(fecundity.m) # random effects summary (between species)

fixefs['lnMass', 'Estimate'] + ranefs$species[, 'Estimate', 'lnMass'] # Slope - species level coeffienct random effects + deviation from mean
exp(fixefs['Intercept', 'Estimate'] + ranefs$animal[, 'Estimate', 'Intercept'] + ranefs$species[, 'Estimate', 'Intercept']) # intercept - estimated on the log scale, take the exponential of the sum to get actual values

# extract the full posterior distribution of species-level intercepts and slopes
posteriors  <-  brms::posterior_samples(fecundity.m)
nrow(posteriors)
dim(posteriors) # 1e4 posterior samples
head(posteriors) 
# b_* indicates the fixed effects; 
# r_animal* intercept-level phylogenetic random effects; 
# r_species* intercept- and slope-level within-species random effects;

# check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
library(bayesplot)
library(tidyverse)
turtle.dat2 <- turtle.dat %>% drop_na(lnClutch, lnMass)
y      <- turtle.dat2$lnClutch
length(y)
yrep   <- posterior_predict(fecundity.m, draws = 500)
dim(yrep)
group  <- turtle.dat2$species
ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)


# for each row in the posterior matrix, obtain the "animal" random effect estimates for the intercept as well as the "species" random effects on the intercept and slope
# intercept = random_intercept_animal + random_intercept_species + fixed_intercept
# slope = random_slope_species + fixed_slope
out_list  <-  vector(mode = 'list', length = length(species))

for (j in seq_along(species)) {
  out_list[[j]]  <-  data.frame(matrix(NA, nrow(posteriors), 2))
  dat_names      <-  paste0(c('Intercept_', 'Slope_'), species[j])
  names(out_list[[j]])  <-  dat_names
  for (i in 1:nrow(posteriors)) {
    out_list[[j]][i, dat_names[1]]  <-  exp(posteriors[i, 'b_Intercept'] + posteriors[i, paste0('r_animal[', species[j], ',Intercept]')] + posteriors[i, paste0('r_species[', species[j], ',Intercept]')])
    out_list[[j]][i, dat_names[2]]  <-  posteriors[i, 'b_lnMass'] + posteriors[i, paste0('r_species[', species[j], ',lnMass]')]
  }
}
out_data  <-  do.call('cbind.data.frame', out_list)
head(out_data)
sapply(out_data, hist) # posterior distributions of coefficients, you can probably make it prettier

# extract 4000 random rows for reproductive output analysis
fecundity.m2 <- sample_n(out_data, 4000)

# Calculate egg.volume slope
turtle.dat$lnVol  <-  log(turtle.dat$egg.volume)

# Test for phylogenetic effects on the intercept, but also a repeated measures effect (i.e. within-species variance) on both the intercept and slope
eggvolume.m <- brm(
  lnVol ~ lnMass + (1 | animal) + (1 + lnMass | species), data = turtle.dat, 
  family = gaussian(), cov_ranef = list(animal = phylo_cor),
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  ),
  chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(eggvolume.m) # show summary of results

plot(eggvolume.m, N = 2, ask = FALSE) # Histograms of posterior samples and trace plots of the intercept
plot(conditional_effects(eggvolume.m), points = TRUE) 

# estimate phylogenetic signal (delta) via hypothesis method
hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(eggvolume.m, hyp, class = NULL))
# delta = 0.78  under estimates

get_variables(eggvolume.m) # list of obtainable variables

fixefs.egg  <-  fixef(eggvolume.m) # fixed effects summary (all data)
ranefs.egg  <-  ranef(eggvolume.m) # random effects summary (between species)

fixefs.egg['lnMass', 'Estimate'] + ranefs.egg$species[, 'Estimate', 'lnMass']
exp(fixefs.egg['Intercept', 'Estimate'] + ranefs.egg$animal[, 'Estimate', 'Intercept'] + ranefs.egg$species[, 'Estimate', 'Intercept'])

# extract the full posterior distribution of species-level intercepts and slopes
posteriors.egg  <-  brms::posterior_samples(eggvolume.m)
dim(posteriors.egg) # 1e4 posterior samples
head(posteriors.egg) 

# check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
turtle.dat2 <- turtle.dat %>% drop_na(lnVol, lnMass)
y      <- turtle.dat2$lnVol
length(y)
yrep   <- posterior_predict(eggvolume.m, draws = 500)
group  <- turtle.dat2$species
ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)

out_list  <-  vector(mode = 'list', length = length(species))

for (j in seq_along(species)) {
  out_list[[j]]  <-  data.frame(matrix(NA, nrow(posteriors.egg), 2))
  dat_names      <-  paste0(c('Intercept_', 'Slope_'), species[j])
  names(out_list[[j]])  <-  dat_names
  for (i in 1:nrow(posteriors.egg)) {
    out_list[[j]][i, dat_names[1]]  <-  exp(posteriors.egg[i, 'b_Intercept'] + posteriors.egg[i, paste0('r_animal[', species[j], ',Intercept]')] + posteriors.egg[i, paste0('r_species[', species[j], ',Intercept]')])
    out_list[[j]][i, dat_names[2]]  <-  posteriors.egg[i, 'b_lnMass'] + posteriors.egg[i, paste0('r_species[', species[j], ',lnMass]')]
  }
}
out_data.egg  <-  do.call('cbind.data.frame', out_list)
head(out_data.egg)
sapply(out_data.egg, hist) # posterior distributions of coefficients, you can probably make it prettier

# extract 4000 random rows for reproductive output analysis
eggvolume.m2 <- sample_n(out_data.egg, 4000)
head(eggvolume.m2)

# combine 4000 random iterations from fecundity.m2 and eggvolume.m2 together
str(fecundity.m2)
library(reshape2)

# merge fecundity intercept data
fecundity.int <- melt(fecundity.m2, measure.vars = c("Intercept_Caretta_caretta","Intercept_Chelonia_mydas",
                                                    "Intercept_Dermochelys_coriacea", "Intercept_Eretmochelys_imbricata",
                                                    "Intercept_Lepidochelys_kempii", "Intercept_Lepidochelys_olivacea",
                                                    "Intercept_Natator_depressus"),
                      variable.name = "species", value.name = "fecundity.int")

# merge fecundity slope data
fecundity.slope <- melt(fecundity.m2, measure.vars = c("Slope_Caretta_caretta","Slope_Chelonia_mydas",
                                                     "Slope_Dermochelys_coriacea", "Slope_Eretmochelys_imbricata",
                                                     "Slope_Lepidochelys_kempii", "Slope_Lepidochelys_olivacea",
                                                     "Slope_Natator_depressus"),
                      variable.name = "species", value.name = "fecundity.slope")

fecundity.sum <- data.frame(fecundity.int$species, fecundity.int$fecundity.int, fecundity.slope$fecundity.slope)
head(fecundity.sum)
names(fecundity.sum)[1] <- "species"
names(fecundity.sum)[2] <- "fecundity.int"
names(fecundity.sum)[3] <- "fecundity.slope"

fecundity.sum$species <- as.factor(substr(fecundity.sum$species, 11, 100)) #remove the word "intercept"
str(fecundity.sum)

# merge egg volume intercept data
eggvolume.int <- melt(eggvolume.m2, measure.vars = c("Intercept_Caretta_caretta","Intercept_Chelonia_mydas",
                                                     "Intercept_Dermochelys_coriacea", "Intercept_Eretmochelys_imbricata",
                                                     "Intercept_Lepidochelys_kempii", "Intercept_Lepidochelys_olivacea",
                                                     "Intercept_Natator_depressus"),
                      variable.name = "species", value.name = "eggvol.int")

# merge egg volume slope data
eggvolume.slope <- melt(eggvolume.m2, measure.vars = c("Slope_Caretta_caretta","Slope_Chelonia_mydas",
                                                      "Slope_Dermochelys_coriacea", "Slope_Eretmochelys_imbricata",
                                                      "Slope_Lepidochelys_kempii", "Slope_Lepidochelys_olivacea",
                                                      "Slope_Natator_depressus"),
                       variable.name = "species", value.name = "eggvol.slope")

eggvolume.sum <- data.frame(eggvolume.int$species, eggvolume.int$eggvol.int, eggvolume.slope$eggvol.slope)
head(eggvolume.sum)
names(eggvolume.sum)[1] <- "species"
names(eggvolume.sum)[2] <- "eggvol.int"
names(eggvolume.sum)[3] <- "eggvol.slope"

eggvolume.sum$species <- as.factor(substr(eggvolume.sum$species, 11, 100)) #remove the word "intercept"
str(eggvolume.sum)

# merge fecundity and egg volume data
slope.m <- data.frame(fecundity.sum, eggvolume.sum) # change to data.frame
str(slope.m)

slope.m$clutchvolume.int <- slope.m$fecundity.int * slope.m$eggvol.int # calculate unlog intercept

slope.m$clutchvolume.slope <- slope.m$fecundity.slope + slope.m$eggvol.slope # calculate clutch volume slope for each species
str(slope.m)

alpha <-  0.05
slope.sum <- slope.m %>% 
  group_by(species) %>% 
  summarise(mean = mean(clutchvolume.slope),
            sd = sd(clutchvolume.slope),
            se = sd(clutchvolume.slope) / sqrt((length(clutchvolume.slope))),
            lower = mean(clutchvolume.slope) - qt(1- alpha/2, (n() - 1))*sd(clutchvolume.slope)/sqrt(n()),
            upper = mean(clutchvolume.slope) + qt(1- alpha/2, (n() - 1))*sd(clutchvolume.slope)/sqrt(n()))

slope.sum <- as.data.frame(slope.sum)
row.names(slope.sum) <- slope.sum$species


# subset to 1000 iterations for graphing
slope.m2 <- sample_n(slope.m, 1000)
str(slope.m2)

# calculate egg energy slope (kJ)
energy.dat <- read.csv("egg energy.csv") # open egg energy file
str(energy.dat)

energy.dat$lnEnergy <- log(energy.dat$egg.energy)
energy.dat$lnMass   <- log(energy.dat$egg.mass)
energy.dat$animal    <-  as.factor(energy.dat$species) # rename animal for phylogenetic random effects, vs within species effect

eggenergy.m <- brm(lnEnergy ~ lnMass, data = energy.dat,
                   chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3
                   )
summary(eggenergy.m)
# fixed effects summary
fixef(eggenergy.m)

# check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
y      <- energy.dat2$lnEnergy
length(y)
yrep   <- posterior_predict(eggenergy.m, draws = 500)
group  <- energy.dat2$species
ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)
ppc_scatter_avg(y, yrep, alpha = 0.7)

# check colour-blind friendly palettes
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
mycol <- brewer.pal(7, "Dark2")
mycol <- c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")


# for subset of 1000 iterations for fecundity
ggplot(turtle.dat, aes(x = lnMass, y = lnClutch, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_abline(data = slope.m2, aes(intercept = log(fecundity.int), slope = fecundity.slope, colour = species), alpha = 0.2) +
  xlab(expression("Female body mass (10"^"x"*", kg)")) + 
  ylab(expression("Fecundity (total clutch output, 10"^"y"*")")) +
  scale_colour_manual(values = mycol)

# for subset of 1000 iterations for egg volume
ggplot(turtle.dat, aes(x = lnMass, y = lnVol, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_abline(data = slope.m2, aes(intercept = log(eggvol.int), slope = eggvol.slope, colour = species), alpha = 0.2) +
  xlab(expression("Female body mass (10"^"x"*", kg)")) + 
  ylab(expression("Egg volume (10"^"y"*", cm"^"3"*")")) +
  scale_colour_manual(values = mycol)


# mean fecundity
fecundity.sum2 <- fecundity.sum %>%
  group_by(species) %>% 
  summarise(mean.int = mean(log(fecundity.int)),
            mean.slope = mean(fecundity.slope))

mass.plot <- ggplot(turtle.dat, aes(x = lnMass, y = lnClutch, fill = species, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_abline(data = fecundity.sum2, aes(intercept = mean.int, slope = mean.slope, colour = species), size = 1.5) +
  xlab(expression("Female body mass (10"^"x"*", kg)")) + 
  ylab(expression("Fecundity (total clutch output, 10"^"y"*")")) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol)

# mean egg volume
eggvolume.sum2 <- eggvolume.sum %>%
  group_by(species) %>% 
  summarise(mean.int = mean(log(eggvol.int)),
            mean.slope = mean(eggvol.slope))

volume.plot <- ggplot(turtle.dat, aes(x = lnMass, y = lnVol, fill = species, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_abline(data = eggvolume.sum2, aes(intercept = mean.int, slope = mean.slope, colour = species), size = 1.5) +
  xlab(expression("Female body mass (10"^"x"*", kg)")) + 
  ylab(expression("Egg volume (10"^"y"*", cm"^"3"*")")) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol)


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

str(slope.m)

tree <- phylo_branch
setdiff(slope.sum$species, as.character(tree$tip.label)) # listed in the database but not in the tree
setdiff(as.character(tree$tip.label), slope.sum$species) # listed in the tree but not in the database
# all ok!

both <- intersect(tree$tip.label, rownames(slope.sum)) # matching tree tip label with rownames
tree2 <- drop.tip(tree, setdiff(both, tree$tip.label)) # drop any tips not matching rownames
data <- slope.sum[both,] # attach intersect with dataframe

slope <- as.matrix(slope.sum)[,2] # extract trait (clutch volume slope)
fit <- phytools::fastAnc(tree2, slope, vars = TRUE, CI = TRUE) # calculate fast estimation of ML ancestral states

td <- data.frame(node = nodeid(tree2, names(slope)),
                 trait = slope) # combine length with tree nodes
td$trait <- as.numeric(levels(td$trait))[td$trait] # change to numeric
nd <- data.frame(node = names(fit$ace), trait = fit$ace) # combine ancestral states with tree nodes

d <- rbind(td, nd)
d$node <- as.numeric(d$node) # change to numeric
tree3 <- full_join(tree2, d, by = 'node') # join trait with tree

# plot tree
tree.plot <- ggtree(tree3, ladderize = FALSE, size = 1) +
  geom_tiplab(size = 3, hjust = -0.1)

# combine with slope points
d1 <- data.frame(species = tree2$tip.label, clutchvolume.slope = data$mean, xmin = data$lower, xmax = data$upper)

slope.plot <- facet_plot(tree.plot + xlim_tree(1.3), panel = "Clutch volume slope",
                         data = d1, geom = geom_pointrangeh, 
                         mapping = aes(x = clutchvolume.slope, xmin = xmin, xmax = xmax, colour = species), size = 1) +
  theme_tree2() + mytheme +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = mycol)


# phylogenetic signal
trait <- slope.sum[,2] # extract clutch volume slope
names(trait) <- slope.sum[,1] # extract species name

# test Pagels lambda with 999 randomizations:
phylosig(phylo_branch, trait, method = "lambda", test = TRUE, nsim = 999)

plot_grid(output.plot, slope.plot, labels = c("", "C"), nrow = 2, rel_heights = c(1, 0.7))

## WITHIN C. MYDAS POPULATION ##-----------------------------------------------------------------------------------
nest.dat <- read.csv("nesting year_2.csv")
str(nest.dat)
names(nest.dat)

length(unique(nest.dat$ID)) # number of individuals tagged
min(nest.dat$sum_n.of.eggs[!is.na(nest.dat$sum_n.of.eggs)])

nest.sum <- nest.dat %>% 
  group_by(year, recruit) %>% 
  summarise(count = length(unique(ID[!is.na(ID)])),
            mass = mean(mass[!is.na(mass)]),
            clutch = mean(sum_n.of.eggs[!is.na(sum_n.of.eggs)])
  )

count(nest.dat, vars = c("ID", "year", "recruit"))

ggplot(nest.sum, aes(x = year, y = count, fill = recruit)) +
  mytheme + geom_line() +
  geom_point(size = 3, shape = 21, colour = "black") + 
  scale_fill_discrete_sequential(palette = "BluGrn")


library(ggpointdensity) # plot density in ggplot2
library(quantreg) # plot quartile regression

# calculate fercundity slope
nest.dat$lnEggs <- log(nest.dat$sum_n.of.eggs)
nest.dat$lnMass <- log(nest.dat$mass)

# calculate quartile range for 
m_quantile_10  <-  brm(bf(lnEggs ~ lnMass, quantile = 0.1), data = nest.dat, family = asym_laplace())
m_quantile_50  <-  brm(bf(lnEggs ~ lnMass, quantile = 0.5), data = nest.dat, family = asym_laplace())
m_quantile_90  <-  brm(bf(lnEggs ~ lnMass, quantile = 0.9), data = nest.dat, family = asym_laplace())
fecund.m       <-  brm(lnEggs ~ lnMass, data = nest.dat)

summary(fecund.m)
nest.dat2 <- nest.dat %>% drop_na(lnEggs)
y      <- nest.dat2$lnEggs
yrep   <- posterior_predict(m_quantile_90, draws = 500)
ppc_scatter_avg(y, yrep, alpha = 0.7)

# fixed effects summary
quant10.m  <-  fixef(m_quantile_10)
quant50.m  <-  fixef(m_quantile_50) 
quant90.m  <-  fixef(m_quantile_90) 
fecund.m2  <-  fixef(fecund.m) 

exp.q10    <-  exp(quant10.m["Intercept","Estimate"])
exp.q50    <-  exp(quant50.m["Intercept","Estimate"])
exp.q90    <-  exp(quant90.m["Intercept","Estimate"])

# plot fecundity with q10-90
nest.plot <- ggplot(nest.dat, aes(x = lnMass, y = lnEggs)) + mytheme +
  geom_pointdensity(adjust = 4, size = 3) +
  geom_abline(intercept = quant10.m["Intercept","Estimate"], 
              slope = quant10.m["lnMass","Estimate"], colour = "#9AD4A8", size = 2) +
  geom_abline(intercept = quant50.m["Intercept","Estimate"], 
              slope = quant50.m["lnMass","Estimate"], colour = "#4CA38F", size = 2) +
  geom_abline(intercept = quant90.m["Intercept","Estimate"], 
              slope = quant90.m["lnMass","Estimate"], colour = "#14505C", size = 2) +
  ylab(expression(paste("Fecundity (total clutch yr"^"-1"*")"))) + 
  xlab("Female mass, Mi, (kg)") +
  scale_colour_continuous_sequential(palette = "BluGrn")

# extract posteriors
quant10.post  <-  brms::posterior_samples(m_quantile_10)
quant50.post  <-  brms::posterior_samples(m_quantile_50)
quant90.post  <-  brms::posterior_samples(m_quantile_90)
fecund.post   <-  brms::posterior_samples(fecund.m)

head(quant90.post)

# calculate egg volume slope - from species dataset

# calculate egg energy slope (kJ)
energy.dat <- read.csv("egg energy.csv") # open egg energy file
str(energy.dat)

energy.dat$lnEnergy <- log(energy.dat$egg.energy)
energy.dat$lnMass   <- log(energy.dat$egg.mass)
energy.dat2 <- subset(energy.dat, species == "Chelonia_mydas")

eggenergy.m <- brm(lnEnergy ~ lnMass, data = energy.dat2)
summary(eggenergy.m)

# fixed effects summary
fixef(eggenergy.m)

energy.post  <-  brms::posterior_samples(eggenergy.m) # egg energy

# extract 4000 random rows for reproductive output analysis
energy.post2 <- sample_n(energy.post, 4000)
head(energy.post2)

# predict expected egg energy for a female of given mass
# Posteriors
str(energy.post2) # egg energy
str(quant10.post) # fecundity 10%
str(quant50.post) # fecundity 50%
str(quant90.post) # fecundity 90%
str(fecund.post)

str(slope.m) 
c.mydas.slope <- subset(slope.m, species %in% "Chelonia_mydas")
str(c.mydas.slope) # fecundity, egg volume, clutch volume

c.mydas.slope$energy.int         <- exp(energy.post2$b_Intercept)
c.mydas.slope$energy.slope       <- energy.post2$b_lnMass
c.mydas.slope$fecundity10.int    <- exp(quant10.post$b_Intercept)
c.mydas.slope$fecundity10.slope  <- quant10.post$b_lnMass
c.mydas.slope$fecundity50.int    <- exp(quant50.post$b_Intercept)
c.mydas.slope$fecundity50.slope  <- quant50.post$b_lnMass
c.mydas.slope$fecundity90.int    <- exp(quant90.post$b_Intercept)
c.mydas.slope$fecundity90.slope  <- quant90.post$b_lnMass
c.mydas.slope$fecund.int         <- exp(fecund.post$b_Intercept)
c.mydas.slope$fecund.slope       <- fecund.post$b_lnMass

# calculate egg energy per volume
c.mydas.slope$eggenergy.int <- c.mydas.slope$eggvol.int * c.mydas.slope$energy.int # calculate log intercept (convert gram to kg)

c.mydas.slope$eggenergy.slope <- c.mydas.slope$eggvol.slope * c.mydas.slope$energy.slope # calculate unlog slope

# calculate total energy output
# 10 %
int10.m <- mean(c.mydas.slope$fecundity10.int * c.mydas.slope$eggenergy.int)
slope10.m <- mean(c.mydas.slope$fecundity10.slope + c.mydas.slope$eggenergy.slope)

# 50 %
int50.m <- mean(c.mydas.slope$fecundity50.int * c.mydas.slope$eggenergy.int)
slope50.m <- mean(c.mydas.slope$fecundity50.slope + c.mydas.slope$eggenergy.slope)

# 90 %
int90.m <- mean(c.mydas.slope$fecundity90.int * c.mydas.slope$eggenergy.int)
slope90.m <- mean(c.mydas.slope$fecundity90.slope + c.mydas.slope$eggenergy.slope)

# mean
int.m <- mean(c.mydas.slope$fecund.int * c.mydas.slope$eggenergy.int)
slope.m <- mean(c.mydas.slope$fecund.slope + c.mydas.slope$eggenergy.slope)

# calculate predicted isometry and allometry relationship for C. mydas
mass <-  seq(min(nest.dat$mass - 10), max(nest.dat$mass + 10)) # in kg
#isometry <- (greenint.m * mass ^ 1 * 0.001) # predicted if isometric, also converted kJ to Mj
hyperallo.10 <- (int10.m * mass ^ slope10.m * 0.001) # from results into MJ
hyperallo.50 <- (int50.m * mass ^ slope50.m * 0.001) # from results into MJ
hyperallo.90 <- (int90.m * mass ^ slope90.m * 0.001) # from results into MJ
hyperallometry <- (int.m * mass ^ slope.m * 0.001) # from results into MJ

# estimate reproductive-energy output based on known mass
int50.m * 76 ^ slope50.m * 0.001

# plot predicted reproductive-energy output
pred.plot <- ggplot() + mytheme +
  #geom_smooth(aes(x = mass, y = isometry), size = 1.5, colour = "black", linetype = "dashed") +
  geom_smooth(aes(x = mass, y = hyperallo.10), size = 1, colour = "#348781", linetype= "dashed") +
  geom_smooth(aes(x = mass, y = hyperallo.50), size = 2, colour = "#348781") +
  #geom_smooth(aes(x = mass, y = hyperallometry), size = 2, colour = "#348781") +
  geom_smooth(aes(x = mass, y = hyperallo.90), size = 1, colour = "#348781", linetype= "dashed") +
  xlab("Female mass, Mi, (kg)") + ylab("Reproductive-energy output, Ri, (MJ)")


## TEMPORAL VARIATION ##--------------------------------------------------------------------------
# changes in mass mass per year
library(plyr) # obtain summary stats
nest.mass <- ddply(nest.dat, "year", summarise,
                   N = length(mass),
                   mean = mean(mass),
                   sd = sd(mass),
                   se = sd/sqrt(N))
names(nest.mass)
nest.mass$output <- int.m * nest.mass$mean ^ slope.m * 0.001

# brm 
yearmass.m <- brm(mass ~ year, data = nest.dat,
  chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3)

summary(yearmass.m)
fixef(yearmass.m)

# posterior check
y      <- nest.dat$mass
yrep   <- posterior_predict(yearmass.m, draws = 500)
ppc_scatter_avg(y, yrep, alpha = 0.7)

# extract posterior samples
yearmass.post  <-  brms::posterior_samples(yearmass.m) # egg energy
head(yearmass.post)

# obtain mean intercept and slope
yearmass.sum <- yearmass.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_year))

# subset 1000 interations
yearmass.post <- sample_n(yearmass.post, 1000) 

yearmass.plot <- ggplot(nest.mass, aes(x = year, y = mean)) + mytheme +
  geom_abline(data = yearmass.post, aes(intercept = b_Intercept, slope = b_year), # 1000 iterations
              alpha = 0.02, colour = "#1B9E77") +
  geom_abline(data = yearmass.sum, aes(intercept = intercept, slope = slope), # mean fit
              lwd = 1, colour = "#1B9E77", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, colour = year), width = 0.2, size = 1) +
  geom_line(linetype = "dashed") +
  geom_point(aes(colour = year), size = 3.5) + xlim(1992, 2020) +
  ylab("Female mass, Mi, (kg)") + xlab("Year") +
  scale_colour_continuous_sequential(palette = "BluGrn")


# n of eggs vs year
nest.egg <- ddply(nest.dat, "year", summarise,
                  N = length(sum_n.of.eggs),
                  mean = mean(sum_n.of.eggs),
                  sd = sd(sum_n.of.eggs),
                  se = sd/sqrt(N))

# brm 
yearegg.m <- brm(sum_n.of.eggs ~ year, data = nest.dat,
  chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3)

summary(yearegg.m)
fixef(yearegg.m)

nest.dat2 <- nest.dat %>% drop_na(sum_n.of.eggs)
y      <- nest.dat2$sum_n.of.eggs
yrep   <- posterior_predict(yearegg.m, draws = 500)
ppc_scatter_avg(y, yrep, alpha = 0.7)

# extract posterior samples
yearegg.post  <-  brms::posterior_samples(yearegg.m) # egg energy
head(yearegg.post)

# obtain mean intercept and slope
yearegg.sum <- yearegg.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_year))

# subset 1000 interations
yearegg.post <- sample_n(yearegg.post, 1000) 

yearegg.plot <- ggplot(nest.egg, aes(x = year, y = mean)) + mytheme +
  geom_abline(data = yearegg.post, aes(intercept = b_Intercept, slope = b_year), # 1000 iterations
              alpha = 0.02, colour = "#1B9E77") +
  geom_abline(data = yearegg.sum, aes(intercept = intercept, slope = slope), # mean fit
              lwd = 1, colour = "#1B9E77", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, colour = year), width = 0.2, size = 1) +
  geom_line(linetype = "dashed") +
  geom_point(aes(colour = year), size = 3.5) + xlim(1992, 2020) +
  ylab("Fecundity") + xlab("Year") +
  scale_colour_continuous_sequential(palette = "BluGrn")

# energy output - compare actual vs predicted
str(nest.dat)
energy <- mean(energy.dat$egg.energy[species == "Chelonia_mydas"], na.rm = TRUE) # mean energy per egg in C. mydas
nest.dat$sumenergy  <- nest.dat$sum_n.of.eggs * energy * 0.001 # mean clutch energy by n of eggs (convert to MJ)
nest.dat$sumenergy2 <- int.m * nest.dat$mass ^ slope.m * 0.001 # predicted mean clutch energy by mass (convert to MJ)
nest.dat$sumenergy10 <- int10.m * nest.dat$mass ^ slope10.m * 0.001 # predicted mean clutch energy by mass (convert to MJ)

nest.dat$lnenergy <- log(nest.dat$sumenergy)
nest.dat$lnenergy2 <- log(nest.dat$sumenergy2)

energy.m <- brm(sumenergy ~ year, data = nest.dat) # actual energy
energy.m2 <- brm(sumenergy2 ~ year, data = nest.dat) # predcited mean energy
energy10.m <- brm(sumenergy10 ~ year, data = nest.dat) # predicted lower quartile

fixef(energy.m)
fixef(energy.m2)
fixef(energy10.m)

# extract posterior samples
energy.post  <-  brms::posterior_samples(energy.m)
energy.post2  <-  brms::posterior_samples(energy.m2)
energy10.post  <-  brms::posterior_samples(energy10.m) 
head(energy.post)

# obtain mean intercept and slope
energy.sum <- energy.post %>%
  summarise(intercept = mean(b_Intercept),
            slope = mean(b_year))

# subset 1000 interations
energy.post <- sample_n(energy.post, 1000) 
energy.post2 <- sample_n(energy.post2, 1000) 
energy10.post <- sample_n(energy10.post, 1000) 

mean(energy10.post$b_Intercept)
mean(energy10.post$b_year)

nest.energy <- ddply(nest.dat, "year", summarise,
                  N = length(sumenergy),
                  mean = mean(sumenergy),
                  sd = sd(sumenergy),
                  se = sd/sqrt(N))

yearenergy.plot <- ggplot(nest.energy, aes(x = year, y = mean)) + mytheme +
  geom_abline(data = energy.post2, aes(intercept = b_Intercept, slope = b_year), # 1000 iterations predict (mean)
              alpha = 0.02, colour = "#858585") +
  geom_abline(data = energy10.post, aes(intercept = b_Intercept, slope = b_year), # 1000 iterations predict (10% quartile)
              alpha = 0.02, colour = "#858585") +
  geom_abline(data = energy.post, aes(intercept = b_Intercept, slope = b_year), # 1000 iterations actual
              alpha = 0.02, colour = "#348781") +
  geom_abline(data = energy.sum, aes(intercept = intercept, slope = slope),
              colour = "#1B9E77") +
  #geom_point(aes(colour = year), size = 3.5) + 
  xlim(1992, 2020) + ylim(50,250) + xlab("Year") +
  ylab("Reproductive-energy output, Ri, (MJ)") +
  scale_colour_continuous_sequential(palette = "BluGrn")

# group plots
plot_grid(pred.plot, yearmass.plot, yearegg.plot, yearenergy.plot, 
          labels = c('A','B','C', 'D'), ncol = 2, align = "v", axis = 'lr')


## SUPPLEMENTARY FILE DATA ##---------------------------------------------------------------------------
# Calculate hatchling slope
head(turtle.dat)
turtle.dat$lnHatch    <-  log(turtle.dat$hatchling.length)

# Test for phylogenetic effects on the intercept, but also a repeated measures effect (i.e. within-species variance) on both the intercept and slope
hatchling.m <- brm(
  lnHatch ~ lnMass + (1 | animal) + (1 + lnMass | species), data = turtle.dat, 
  family = gaussian(), cov_ranef = list(animal = phylo_cor),
  prior = c(
    prior(normal(0, 10), "b"), # mean of 0 and SD of 10 (wide distribution)
    prior(normal(0, 10), "Intercept"), # mean of 0 and SD of 10 (if log then wide negative disribution becomes infinite)
    prior(student_t(3, 0, 20), "sd"), # class of random effect deviation to calcuate - has to be postiive (no negative SD)
    prior(student_t(3, 0, 20), "sigma") # residual SD parameter
  ),
  chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15)
)

summary(hatchling.m) # show summary of results
fixef(hatchling.m)

plot(conditional_effects(hatchling.m), points = TRUE) 

# estimate phylogenetic signal (delta) via hypothesis method
hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(hatchling.m, hyp, class = NULL))
# delta = 0.81 under estimates

# extract the full posterior distribution of species-level intercepts and slopes
posteriors  <-  brms::posterior_samples(hatchling.m)
nrow(posteriors)
dim(posteriors) # 1e4 posterior samples
head(posteriors) 

# check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
turtle.dat2 <- turtle.dat %>% drop_na(lnHatch, lnMass)
y      <- turtle.dat2$lnHatch
yrep   <- posterior_predict(hatchling.m, draws = 500)
group  <- turtle.dat2$species
ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)

out_list  <-  vector(mode = 'list', length = length(species))

for (j in seq_along(species)) {
  out_list[[j]]  <-  data.frame(matrix(NA, nrow(posteriors), 2))
  dat_names      <-  paste0(c('Intercept_', 'Slope_'), species[j])
  names(out_list[[j]])  <-  dat_names
  for (i in 1:nrow(posteriors)) {
    out_list[[j]][i, dat_names[1]]  <-  exp(posteriors[i, 'b_Intercept'] + posteriors[i, paste0('r_animal[', species[j], ',Intercept]')] + posteriors[i, paste0('r_species[', species[j], ',Intercept]')])
    out_list[[j]][i, dat_names[2]]  <-  posteriors[i, 'b_lnMass'] + posteriors[i, paste0('r_species[', species[j], ',lnMass]')]
  }
}
out_data  <-  do.call('cbind.data.frame', out_list)
head(out_data)
sapply(out_data, hist) # posterior distributions of coefficients, you can probably make it prettier

# extract 4000 random rows for reproductive output analysis
hatchling.m2 <- sample_n(out_data, 4000)

# merge fecundity intercept data
hatchling.int <- melt(hatchling.m2, measure.vars = c("Intercept_Caretta_caretta","Intercept_Chelonia_mydas",
                                                     "Intercept_Dermochelys_coriacea", "Intercept_Eretmochelys_imbricata",
                                                     "Intercept_Lepidochelys_kempii", "Intercept_Lepidochelys_olivacea",
                                                     "Intercept_Natator_depressus"),
                      variable.name = "species", value.name = "hatchling.int")

# merge fecundity slope data
hatchling.slope <- melt(hatchling.m2, measure.vars = c("Slope_Caretta_caretta","Slope_Chelonia_mydas",
                                                       "Slope_Dermochelys_coriacea", "Slope_Eretmochelys_imbricata",
                                                       "Slope_Lepidochelys_kempii", "Slope_Lepidochelys_olivacea",
                                                       "Slope_Natator_depressus"),
                        variable.name = "species", value.name = "hatchling.slope")

hatchling.sum <- data.frame(hatchling.int$species, hatchling.int$hatchling.int, hatchling.slope$hatchling.slope)
head(hatchling.sum)
names(hatchling.sum)[1] <- "species"
names(hatchling.sum)[2] <- "hatchling.int"
names(hatchling.sum)[3] <- "hatchling.slope"
hatchling.sum$species <- as.factor(substr(hatchling.sum$species, 11, 100)) #remove the word "intercept"

hatchling.sum2 <- hatchling.sum %>%
  group_by(species) %>% 
  summarise(mean.int = mean(log(hatchling.int)),
            mean.slope = mean(hatchling.slope))

hatch.plot <- ggplot(turtle.dat, aes(x = lnMass, y = lnHatch, fill = species, colour = species)) +
  mytheme + geom_point(size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_abline(data = hatchling.sum2, aes(intercept = mean.int, slope = mean.slope, colour = species), size = 1.5) +
  xlab(expression("Female body mass (10"^"x"*", kg)")) + 
  ylab(expression("Hatchling size (10"^"y"*", mm)")) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol)

# plot changes in body mass per year (distribution plot)
install.packages("ggridges")
library(ggridges)

nest.dat2 <- transform(nest.dat, year2 = as.factor(year)) # for geom_density to work

ggplot(nest.dat2, aes(x = mass, y = year2, fill = recruit)) + mytheme +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.5) +
  xlab("Female mass (kg)") + ylab("Year") +
  scale_fill_discrete_sequential(palette = "BluGrn")




## Extra code ##---------------------------------------------------------------------------------------------
# create sequence of mass values for each species
# Make sure you have a vector with each of your species names
library(data.table)
myids <- unique(turtle.dat$species)

mass.group <- list() # create an empty list to store your results into

for(i in 1:length(myids)){ # for each possible position for the length of the ID list
  # First subset out the group of interest
  IDgroup <- subset(turtle.dat, species == myids[i])
  # for each position i, create a dataframe in the list
  mass.group[[i]] <- data.frame(seq(min(IDgroup$pred.mass - 5, na.rm = TRUE), max(IDgroup$pred.mass + 5, na.rm = TRUE)))
  # Add an extra column to the dataframe you just made with the ID so you know which sequence it was
  mass.group[[i]]$species <- myids[i]
}
# Combine each of the lists you made here into a single dataframe/table
mass.group <- data.table(Reduce(rbind, mass.group))

str(mass.group)
names(mass.group)[1] <- "mass" # rename first column


# for c mydas calculation
# old code from lm (not brm)
eggenergy.m2$intercept <- exp(coef(eggvolume.m)['(Intercept)']) * 
  exp(coef(eggenergy.m)['(Intercept)']) # calculate log intercept

eggenergy.m2$slope <- coef(eggvolume.m)['log(mass)'] * 
  coef(eggenergy.m)['log(egg.mass)'] # calculate unlog slope

# calculate total energy output
greenint.m <- exp(coef(total.nest.m)['(Intercept)']) * eggenergy.m2$intercept
#greenint.m2 <- 10^greenint.m # convert to power function
greenslope.m <- coef(total.nest.m)['log(mass)'] + eggenergy.m2$slope

# how much energy per egg (kJ)
nest.egg$sumenergy <- nest.egg$mean * energy * 0.001 # predicted mean clutch energy per year (convert to MJ)

nest.year <- merge(nest.mass, nest.egg, by = "year") # merge both summary data

# plot predicted and actual reproductive output
yearenergy.plot <- ggplot(nest.year, aes(x = year)) + mytheme +
  geom_line(aes(y = output), colour = "#858585",linetype = "dashed") + # prediction
  geom_line(aes(y = sumenergy), size = 1, colour = "#348781") + # actual energy
  geom_point(aes(y = output), colour = "#858585", shape = 17, size = 3.5) +
  geom_point(aes(y = sumenergy, colour = year), size = 3.5) +
  xlim(1992, 2020) + scale_colour_continuous_sequential(palette = "BluGrn") +
  ylab("Mean reproductive-energy output, Ri, (MJ)") + xlab("Year")

