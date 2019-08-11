install.packages("cowplot")
library(ggplot2)
library(cowplot)
mytheme <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                              axis.text = element_text(size = 11, colour = "black"))

## load file ###################################################################################
nesting<-read.csv("C:/Users/nicho/Dropbox (Personal)/Turtle project/nesting year.csv") # load data
head(nesting) # check column headings

# length vs width graph
plot1A <- ggplot(nesting, aes(x=mean_length, y=mean_width)) + mytheme +
  geom_point(shape=21, size=2, fill="white") + ylim(60,120) + xlim(60,130) +
  ylab("shell width (cm)") + xlab("shell length (cm)")


# length vs n of eggs graph
ggplot(nesting, aes(x=log(mean_length), y=log(mean_n.of.eggs))) + mytheme +
  geom_point(shape=21, size=2, fill="white") +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ylab("log number of eggs laid per female") + xlab(" log body mass (kg)")

# mass vs number of eggs graph
library(scales) # to access break formatting functions e.g. log transformed axis
library(ggpmisc) # use stat_poly_eq function
plot1B <- ggplot(nesting, aes(x=log(mass), y=log(mean_n.of.eggs))) + mytheme +
  geom_point(shape=21, size=2, fill="white") +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ylab("log number of eggs laid per female") + xlab(" log body length (cm)")

ylim(0,200) + xlim(50,200) #extra if not log/log

# check model with equation in graph - for mass
lmodel <- lm(mean_n.of.eggs ~ mass, data = nesting)
summary(lmodel)
lmodel$coefficients

# using quantreg package
install.packages("quantreg")
library(quantreg)
# test for mass data
plot(nesting$mass,nesting$mean_n.of.eggs,cex=.25,type="n")
points(nesting$mass,nesting$mean_n.of.eggs,cex=.5,col="blue")
abline(rq(nesting$mean_n.of.eggs~nesting$mass,tau=.5),col="blue") # 0.5 quartile regression
abline(lm(nesting$mean_n.of.eggs~nesting$mass),lty=2,col="red") # the dreaded ols line
taus <- c(.15,.35,.55,.75,.95)
for( i in 1:length(taus)){
  abline(rq(nesting$mean_n.of.eggs~nesting$mass,tau=taus[i]),col="gray")
}

# log transformed mass data of different quartile ranges
plot(nesting$mass,nesting$mean_n.of.eggs,log="xy", xlim=c(50, 250), ylim=c(20, 200))
taus <- c(0.15,0.35,0.55,0.75,0.95) # set quartile ranges
abline(rq(log10(nesting$mean_n.of.eggs)~log10(nesting$mass),tau=.5),col="blue")
abline(lm(log10(nesting$mean_n.of.eggs)~log10(nesting$mass)),lty = 3,col="red")
for( i in 1:length(taus)){
  abline(rq(log10(nesting$mean_n.of.eggs)~log10(nesting$mass),tau=taus[i]),col="gray")
  }

fit2 <- rq(log10(nesting$mean_n.of.eggs)~log10(nesting$mass),tau=c(.15,.35,.55,.75,.95))
summary(fit2)
fit2$coefficients

# test for length data
plot(nesting$mean_length,nesting$mean_n.of.eggs, main = "length vs n of eggs")
taus <- c(0.15,0.35,0.55,0.75,0.95,0.99)
abline(rq(nesting$mean_n.of.eggs~nesting$mean_length,tau=.5),col="blue")
abline(lm(nesting$mean_n.of.eggs~nesting$mean_length),lty = 3,col="red")
for( i in 1:length(taus)){
  abline(rq(nesting$mean_n.of.eggs~nesting$mean_length,tau=taus[i]),col="gray")
}

fit2a <- rq(nesting$mean_n.of.eggs~nesting$mean_length,tau=c(.15,.35,.55,.75,.95,.99))# model summary for all quartiles
summary(fit2a)
fit2a$coefficients

# log transformed length data
plot(nesting$mean_length,nesting$mean_n.of.eggs,log="xy", main = " log length vs n of eggs")
taus <- c(0.15,0.35,0.55,0.75,0.95,0.99)
abline(rq(log10(nesting$mean_n.of.eggs)~log10(nesting$mean_length),tau=.5),col="blue")
abline(lm(log10(nesting$mean_n.of.eggs)~log10(nesting$mean_length)),lty = 3,col="red")
for( i in 1:length(taus)){
  abline(rq(log10(nesting$mean_n.of.eggs)~log10(nesting$mean_length),tau=taus[i]),col="gray")
}

fit2b <- rq(log10(nesting$mean_n.of.eggs)~log10(nesting$mean_length),tau=c(.15,.35,.55,.75,.95,.99))
summary(fit2b)
fit2b$coefficients

# mass vs year graph
nesting$year <- as.factor(nesting$year) # chnage year to factor

library(plyr) # obtain summary stats 
mass <- ddply(nesting, "year", summarise,
               N = length(mass),
               mean = mean(mass),
               sd = sd(mass),
               se = sd/sqrt(N))

plot1C <- ggplot(length, aes(x=year, y=mean)) + mytheme +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  geom_line(linetype = "dashed") +
  geom_smooth(method=lm, color="black") +
  geom_point(shape=21, size=3, fill="white") + ylim(95,102) + xlim(1992,2020) +
  ylab("mean mass of nesting females (kg)") + xlab("year")

# n of eggs vs year
egg <- ddply(nesting, "year", summarise,
                N = length(mean_n.of.eggs),
                mean = mean(mean_n.of.eggs),
                sd = sd(mean_n.of.eggs),
                se = sd/sqrt(N))
 
plot1D <- ggplot(egg, aes(x=year, y=mean)) + mytheme +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  geom_line(linetype = "dashed") +
  geom_smooth(method=lm, color="black") +
  geom_point(shape=21, size=3, fill="white") + ylim(70,110) + xlim(1992,2020) +
  ylab("mean number of eggs laid") + xlab("year")

# group plots
plot_grid(plot1A, plot1B, plot1C, plot1D, labels = c('A','B','C','D'))


## comparative species ######################################################
repro<-read.csv("C:/Users/nicho/Dropbox (Personal)/Turtle project/reproductive output.csv") # load data
head(repro) # check column headings

# length vs clutch
graph1 <- ggplot(repro, aes(x=female.length, y=clutch.size, colour = species)) + mytheme +
  geom_point(size = 2) +
  geom_smooth(method=lm)

# length vs eggs size
# Remove rows that are missing values for egg.diameter
graph2 <- ggplot(repro, aes(x=female.length, y=egg.diameter, colour = species)) + mytheme +
  geom_point(size = 2) +
  geom_smooth(method=lm)

graph3 <- ggplot(repro, aes(x=female.length, y=hatchling.length, colour = species)) + mytheme +
  geom_point(size = 2) +
  geom_smooth(method=lm)

prow <- plot_grid(
  graph1 + theme(legend.position="none"),
  graph2 + theme(legend.position="none"),
  graph3 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)

legend <- get_legend(graph1 + theme(legend.box.margin = margin(0, 0, 0, 12))) # create some space to the left of the legend

# add the legend to the row we made earlier. Give it one-third of the width of one plot (via rel_widths).
plot_grid(prow, legend, rel_widths = c(3, .4))
