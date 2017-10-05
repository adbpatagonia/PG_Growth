

## load libraries ----
library(ggplot2)
library(lubridate)
library(dplyr)
library(betareg)
library(MuMIn)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')

source("functions/HighstatLibV7.R")
source('functions/panelcor.r')
source("functions/MCMCSupportHighstatV2.R")

age <- read.csv('data/PG_age.csv', header = T, as.is = T)
growth <- read.csv('data/PG_growth.csv', header = T, as.is = T)
maturity <-  read.csv('data/PG_maturity.csv', header = T, as.is = T)

fwd <- 2.9
fwj <- 6.3
fwf <- 8.3



## fix records ----
# One seal is coded as age 3, but everything points to it being much older, therefore remove age
age[which(age$idsex == '20082059F'),'cohortage'] <- NA
# one seal has length = 2300, turn it into NA
growth[which(growth$length == 2300),'length'] <- NA
# at least 3 seals have the wrong weight, turn them to NA
ids <-  c("19940928F", "20042978F", "20062138F")
growth[which(growth$idsex %in% ids), 'weight'] <- NA


## merge data sets ----
data <- merge(growth, maturity, by = 'idsex', all = T)
# get rid of duplicate id sex
data <- data[!duplicated(data$idsex),]


## remove fetal weight ----
data$weight <- ifelse(data$fmaturity %in% 2:3 & data$month == 12, data$weight - fwd, data$weight)
data$weight <- ifelse(data$fmaturity %in% 2:3 & data$month == 1, data$weight - fwj, data$weight)
data$weight <- ifelse(data$fmaturity %in% 2:3 & data$month == 2, data$weight - fwf, data$weight)


library(stringr)

## subset data ----
m <- c(1,2,12)
mat <- 2:9

adf <- na.omit(data[which(data$fmaturity %in% mat & data$sex =='F' & data$month %in% m & data$length > 130), c('idsex', 'month', 'length', 'weight', 'fmaturity')])

adf$ll <- log(adf$length); adf$lw <- log(adf$weight)

adfd <- adf[which(adf$month == 12),]
adfj <- adf[which(adf$month == 1),]
adff <- adf[which(adf$month == 2),]
# adfj <- adf %>%
#   filter(month == 1)
# 
# adff <- adf %>%
#   filter(month == 2)


## fit models and get predicted weight ----
fitd <- lm( lw ~ ll, data = adfd)
adfd$predw <- exp(predict(fitd, type = 'response'))
fitj <- lm( lw ~ ll, data = adfj)
adfj$predw <- exp(predict(fitj, type = 'response'))
fitf <- lm( lw ~ ll, data = adff)
adff$predw <- exp(predict(fitf, type = 'response'))

## bind model outputs ----
adf <- rbind(adfd, adfj, adff)

## obtain condition indices ----
adf$cohyear <- as.numeric(str_sub(adf$idsex, 1, 4))
adf$relcond <- with(adf, weight/predw)

## visualize l-w relation----
pdataj <- expand.grid(month = 1, ll = seq(min(adfj$ll),max(adfj$ll), length.out = 100))
pdatad <- expand.grid(month = 12, ll = seq(min(adfd$ll),max(adfd$ll), length.out = 100))
pdataf <- expand.grid(month = 2, ll = seq(min(adff$ll),max(adff$ll), length.out = 100))
pdatad$lpw <- predict(fitd, newdata = pdatad, type = 'response')  
pdataj$lpw <- predict(fitj, newdata = pdataj, type = 'response')  
pdataf$lpw <- predict(fitf, newdata = pdataf, type = 'response') 

pdata <- rbind(pdataj, pdataf, pdatad)

w <- ggplot(adf, aes(exp(ll), exp(lw), colour = as.factor(month)))
w <- w + geom_point()
w <- w + geom_line(data= pdata, aes(exp(ll), exp(lpw), colour = as.factor(month)))
w <- w + theme_bw()
w

## visualize condition ----
pregn <- 2:3
notpreg <- c(4:7,9)
ep <- 8
adf$preg <- ifelse(adf$fmaturity %in% pregn, 1, 
                   ifelse(adf$fmaturity %in% notpreg, 2,
                          ifelse(adf$fmaturity %in% ep, 3, 0)))

meanc <- adf %>%
  group_by(preg, cohyear) %>%
  summarize(meanc = mean(relcond))

cp <- ggplot(adf, aes( factor(cohyear), relcond))
cp <- cp + geom_violin(trim = FALSE)
cp <- cp + geom_jitter(position = position_jitter(0.1)) + stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")
cp <- cp + xlab("Year") + ylab('Condition')
cp <- cp + theme_bw()
cp <- cp + geom_hline(yintercept = 1)
cp <- cp + facet_grid(preg ~ .)
print(cp)


## fecundity ----
# use same data files as used in Stenson et al paper
fecun <- read.csv('data/fecundity_2014-feb20.csv', header = T)
bootci <- read.csv('data/bootCI_50000.csv',header=T)
names(bootci) <- c('cohort_year', 'feclb','fecub','ablb','abub')

bootci[which(bootci$cohort_year == 1979), 'abub'] <- 0

## calculate biological rates ----
fecun$fecrate <- fecun$preg/fecun$mature
fecun$totpreg <- with(fecun,EP + preg)
fecun$abrate <- fecun$EP/fecun$totpreg
fecun$fecratet1 <- c(NA, fecun$fecrate[1:(length(fecun$fecrate) - 1)])
fecun$abratet1 <- c(NA, fecun$abrate[1:(length(fecun$abrate) - 1)])

## calculate mean relative condition by cohort year ----
meancond <- adf %>%
  group_by(cohyear) %>%
  summarize(meanc = mean(relcond))

names(meancond) <- c('cohort_year', 'meancond')
meancond$cohort_year <- as.integer(meancond$cohort_year)

## bootstrap condition data by cohort year
bootcond <- data.frame(cohort_year = sort(unique(adf$cohyear)), lcc = rep(NA, length(unique(adf$cohyear))), ucc = rep(NA, length(unique(adf$cohyear))))

for(i in sort(unique(adf$cohyear))) {
  dd <- quantile(base::sample(adf$relcond[which(adf$cohyear == i)], size = 100000, replace = T), c(0.025, 0.975))
  bootcond[which(bootcond$cohort_year == i), 'lcc'] <- dd[1]
  bootcond[which(bootcond$cohort_year == i), 'ucc'] <- dd[2]
}
rm(i)

## merge fecun and condition data ----
fecun <- merge(fecun, meancond, by = 'cohort_year', all.x = T)
fecun <- merge(fecun, bootcond, by = 'cohort_year', all.x = T)
fecun <- merge(fecun, bootci, by = 'cohort_year', all.x = T)

## plot fecundity ----
ay <- ggplot(fecun, aes(cohort_year, abrate))
ay <- ay + geom_point()
ay <- ay + geom_pointrange(aes(ymin = ablb, ymax = abub))
ay <- ay + theme_bw()
ay


cy <- ggplot(fecun, aes(cohort_year, meancond))
cy <- cy + geom_point()
cy <- cy + geom_pointrange(aes(ymin = lcc, ymax = ucc))
#f <- f + geom_line(aes(cohort_year, abrate))
cy <- cy + theme_bw()
cy

fy <- ggplot(fecun, aes(cohort_year, fecrate))
fy <- fy + geom_point()
fy <- fy + geom_pointrange(aes(ymin = feclb, ymax = fecub))
fy <- fy + theme_bw()
fy


1 - log(fecun$meancond)




fit <- nls(abrate ~ 1-b1*log(meancond), data=fecun, start=list(b1=1))
xdat <- data.frame(meandcond = seq(min(na.omit(fecun$meancond)), max(na.omit(fecun$meancond)), length.out = 100))
xdat$Eabrate <- predict(fit, newdata = xdat)

plot(log(abrate) ~ (meancond), data = fecun, pch =16)

ff <- na.omit(fecun[,c('abrate','meancond'),])
## check 1986!! ----
## ff[which(ff$abrate == 0 & ff$meancond < 0.9), 'meancond'] <- NA
ff$abrate[which(ff$abrate == 0)] <- 1e-16
exponential.model <- lm(log(abrate) ~  meancond, data = ff)
condvalues <- seq(0, 1.5, length.out = 100)
abrate.exponential <- exp(predict(exponential.model,list(meancond=condvalues)))
plot(ff$meancond, ff$abrate, pch=16, xlim = c(0.8, 1.2))
lines(condvalues, abrate.exponential,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")





#### start betaregs ----
fecun$Capt1<-c(NA,fecun$Cap[1:(length(fecun$Cap)-1)])
fecun$Capt2<-c(NA,NA,fecun$Cap[1:(length(fecun$Cap)-2)])
fecun$ArcCodt1<-c(NA,fecun$ArcCod[1:(length(fecun$ArcCod)-1)])
fecun$ArcCodt2<-c(NA,NA,fecun$ArcCod[1:(length(fecun$ArcCod)-2)])
fecun$Sandt1<-c(NA,fecun$Sand[1:(length(fecun$Sand)-1)])
fecun$Sandt2<-c(NA,NA,fecun$Sand[1:(length(fecun$Sand)-2)])
fecunice <- fecun
fecunice <- subset(fecunice, cohort_year<2014 & cohort_year>1995)
fecunice <- fecunice[order(fecunice$cohort_year),]
fecunice$abrate<-replace(fecunice$abrate,fecunice$abrate==0,1e-6)




pairs(c(fecunice[,c('ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1', 'meancond')]),upper.panel=points, lower.panel=panelcor,diag.panel=panel.hist)
fecunicet1 <- fecunice[,c('abrate','ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','cohort_year', 'meancond')]

#Data exploration
#Outliers
MyVar <- c('abrate','ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','meancond')
MyXVar <- c('ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','meancond')
Mydotplot(fecunicet1[,MyVar])

#Collinearity
pairs(fecunicet1[,MyVar], lower.panel = panel.cor)

corvif(fecunicet1[ ,MyXVar])


VIF <-  corvif(fecunicet1[ ,MyXVar])
#setwd(outpath)
#write.csv(VIF,'VIF_abrate_review2.csv', row.names=T) 
#setwd(inpath)




##
### This code constructs and fits all submodels starting from a full model - main effects only 
##
## define link function
linkfn <- 'loglog'
## obtain column names
Cols <- names(fecunicet1)
## exclude response variable, weights and year
Cols <- Cols[! Cols %in% c('cohort_year','weight','abrate')]
n <- length(Cols)
## obtain all combinations of explanatory variables
id <- unlist(
  lapply(1:n,
         function(i)combn(1:n,i,simplify=F)
  )
  ,recursive=F)

## apply these combinations into a formula, with abrate as response variable
Formulas <- sapply(id,function(i)
  paste("abrate~",paste(Cols[i],collapse="+"))
)
## number of models to be fitted            
nmodels <-  length(Formulas)
## fit all models
allmodelst1 <- lapply(Formulas,function(i)
  betareg(as.formula(i),data=fecunicet1,link=linkfn))

## obtain statistics for all models, including the expression of each model    
mst1 <- data.frame(modelnum=as.numeric(1:nmodels))
for (i in 1:nrow(mst1)){
  mst1$model[i] <- as.character(allmodelst1[[i]]$formula[3])
  mst1$N[i] <-   allmodelst1[[i]]$nobs
  mst1$K[i] <- length(allmodelst1[[i]]$coefficients$mean)
  mst1$pseudorsquared[i] <- round(allmodelst1[[i]]$pseudo.r.squared,4)
  mst1$AICc[i] <- AICc(allmodelst1[[i]])
  mst1$LH[i] <- (allmodelst1[[i]]$loglik)
}

### compute AICc and other stats for all models 
mst1$deltaAICc <- mst1$AICc - min(mst1$AICc)

mst1$wi<-exp(-0.5*mst1$deltaAICc)/sum(exp(-0.5*mst1$deltaAICc))
mst1$er <- max(mst1$wi)/mst1$wi
mst1<-mst1[order(mst1$modelnum),]
mst1<-mst1[order(mst1$AICc),]  

mst1[which(mst1$er<10),c('model','deltaAICc','er')]


#Cols <- c(names(fecunicet1tf), names(fecunicet1))
Cols <- c( names(fecunicet1))
Cols <- unique(Cols[! Cols %in% c('cohort_year','weight','abrate')])
relimp <- data.frame(variable=as.character(rep(NA,length(Cols))), wij=as.numeric(rep(NA,length(Cols))), stringsAsFactors = F)
for (i in 1:length(Cols)){
  relimp$variable[i] <- Cols[i]
  relimp$wij[i] <- sum(mst1[grep(Cols[i],mst1$model),'wi'])
}
relimp <- relimp[rev(order(relimp$wij)),]


