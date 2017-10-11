
## load libraries ----
#library(ggplot2)
library(lubridate)
library(dplyr)
library(cowplot)
library(MuMIn)
library(betareg)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')
source("functions/HighstatLibV7.R")
source('functions/panelcor.r')
source("functions/MCMCSupportHighstatV2.R")

load('Rdata/adf.Rdata')

## link function ----
linkfn <- 'log'

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
ncond <- adf %>%
  group_by(cohyear) %>%
  group_size()
meancond <- cbind(meancond, ncond)
rm(ncond)
meancond[which(meancond$ncond < 20), 'meancond'] <- NA


## bootstrap condition data by cohort year
bootcond <- data.frame(cohort_year = sort(unique(adf$cohyear)), lcc = rep(NA, length(unique(adf$cohyear))), ucc = rep(NA, length(unique(adf$cohyear))))
for (i in sort(unique(adf$cohyear))) {
  dd <- quantile(base::sample(adf$relcond[which(adf$cohyear == i)], size = 100000, replace = T), c(0.025, 0.975))
  bootcond[which(bootcond$cohort_year == i), 'lcc'] <- dd[1]
  bootcond[which(bootcond$cohort_year == i), 'ucc'] <- dd[2]
}
rm(i)

## merge fecun and condition data ----
fecun <- merge(fecun, meancond, by = 'cohort_year', all.x = T)
fecun <- merge(fecun, bootcond, by = 'cohort_year', all.x = T)
fecun <- merge(fecun, bootci, by = 'cohort_year', all.x = T)


#rm(bootcond, meancond, bootci, adf)

## plot fecundity ----
ay <- ggplot(fecun, aes(cohort_year, abrate))
ay <- ay + geom_point()
ay <- ay + geom_linerange(aes(ymin = ablb, ymax = abub))
ay <- ay + theme_set(theme_cowplot())
ay <- ay + labs(x = 'Cohort year', y = 'Abortion rate')
ay


cy <- ggplot(fecun, aes(cohort_year, meancond))
cy <- cy + geom_point()
cy <- cy + geom_linerange(data = na.omit(fecun[,c('cohort_year', 'meancond', 'lcc', 'ucc'),]), aes(ymin = lcc, ymax = ucc))
#f <- f + geom_line(aes(cohort_year, abrate))
cy <- cy + theme_set(theme_cowplot())
cy <- cy + labs(x = 'Cohort year', y = 'Relative condition')
cy

fy <- ggplot(fecun, aes(cohort_year, fecrate))
fy <- fy + geom_point()
fy <- fy + geom_linerange(aes(ymin = feclb, ymax = fecub))
fy <- fy + theme_set(theme_cowplot())
fy <- fy + labs(x = 'Cohort year', y = 'Pregnancy rate')
fy

ca <- ggplot(fecun[which(fecun$cohort_year > 1986),], aes(meancond, abrate))
ca <- ca + geom_point() 
ca <- ca + geom_errorbarh(data = fecun, aes(xmin = lcc, xmax = ucc), height = 0)
ca <- ca + geom_errorbar(data = fecun, aes(ymin = ablb, ymax = abub), width = 0)
ca <- ca + ylim(0, 0.6)
ca <- ca + xlim(0.6, 1.6)
ca <- ca + theme_set(theme_cowplot())
ca <- ca + labs(x = 'Relative condition', y = 'Abortion rate')
ca

ca <- ggplot(fecun, aes(meancond, abrate))
ca <- ca + geom_point() 
#ca <- ca + geom_errorbarh(data = fecun, aes(xmin = lcc, xmax = ucc), height = 0)
ca <- ca + geom_errorbar(data = fecun, aes(ymin = ablb, ymax = abub), width = 0)
ca <- ca + ylim(0, 0.6)
#ca <- ca + xlim(0.8, 1.2)
ca <- ca + theme_set(theme_cowplot())
ca <- ca + labs(x = 'Relative condition', y = 'Abortion rate')
ca

# fit <- nls(abrate ~ 1-b1*log(meancond), data=fecun, start=list(b1=1))
# xdat <- data.frame(meandcond = seq(min(na.omit(fecun$meancond)), max(na.omit(fecun$meancond)), length.out = 100))
# xdat$Eabrate <- predict(fit, newdata = xdat)


ff <- na.omit(fecun[,c('abrate','meancond'),])
## check 1986!! ----
## ff[which(ff$abrate == 0 & ff$meancond < 0.9), 'meancond'] <- NA
ff$abrate[which(ff$abrate == 0)] <- 1e-16
exponential.model <- lm(log(abrate) ~  meancond, data = ff)

exponential.models <- lm(abrate ~  exp(1-meancond), data = ff)

condvalues <- data.frame(condvalues = seq(0, 1.5, length.out = 100))
condvalues$Eabrate <- exp(predict(exponential.model,list(meancond = condvalues$condvalues)))
condvalues$Eabrates <- (predict(exponential.models,list(meancond = condvalues$condvalues)))
plot(ff$meancond, ff$abrate, pch = 16, xlim = c(0.8, 1.2))
lines(condvalues$condvalues, condvalues$Eabrate,lwd = 2, col = "red")
lines(condvalues$condvalues, condvalues$Eabrates,lwd = 2, col = "red")

fecun$Eabrate <- exp(predict(exponential.model,list(meancond = fecun$meancond)))





#### start betaregs ----
fecun$Capt1 <- c(NA, fecun$Cap[1:(length(fecun$Cap) - 1)])
fecun$Capt2 <- c(NA, NA, fecun$Cap[1:(length(fecun$Cap) - 2)])
fecun$ArcCodt1 <- c(NA, fecun$ArcCod[1:(length(fecun$ArcCod) - 1)])
fecun$ArcCodt2 <- c(NA, NA, fecun$ArcCod[1:(length(fecun$ArcCod) - 2)])
fecun$Sandt1 <- c(NA, fecun$Sand[1:(length(fecun$Sand) - 1)])
fecun$Sandt2 <- c(NA, NA, fecun$Sand[1:(length(fecun$Sand) - 2)])



fecunice <- fecun
#fecunice <- subset(fecunice, cohort_year < 2014 )

fecunice <- fecunice[,c('cohort_year','fecrate','totpop','abrate','meancond')]
fecunice <- na.omit(fecunice)
########################################################
#Data exploration
#Outliers
MyVar <- c('cohort_year','fecrate','totpop','abrate','meancond')
MyXVar <- c('totpop','abrate','meancond')
#Mydotplot(fecunice[,MyVar])

#Collinearity
pairs(fecunice[,MyVar], lower.panel = panel.cor)

corvif(fecunice[ ,MyXVar])

VIF <-  corvif(fecunice[ ,MyXVar])
#write.csv(VIF,'VIF_fecrate1.csv', row.names=T)
pairs(fecunice[ ,MyXVar], lower.panel = panel.cor)




### This code constructs and fits all submodels starting from a full model - main effects only ----
## obtain column names
Cols <- names(fecunice)
## exclude response variable, weights and year
Cols <- Cols[!Cols %in% c('cohort_year', 'fecrate')]
n <- length(Cols)
## obtain all combinations of explanatory variables
id <- unlist(
  lapply(1:n,
         function(i)combn(1:n,i,simplify = F)
  )
  ,recursive = F)

## apply these combinations into a formula, with fecrate as response variable
Formulas <- sapply(id,function(i)
  paste("fecrate~",paste(Cols[i],collapse = "+"))
)
## number of models to be fitted            
nmodels <-  length(Formulas)
## fit all models ----
allmodels <- lapply(Formulas,function(i)
  betareg(as.formula(i),data = fecunice,link = linkfn))

## obtain statistics for all models, including the expression of each model   ---- 
ms <- data.frame(modelnum = as.numeric(1:nmodels))
for (i in 1:nrow(ms)) {
  ms$model[i] <- as.character(allmodels[[i]]$formula[3])
  ms$N[i] <-   allmodels[[i]]$nobs
  ms$K[i] <- length(allmodels[[i]]$coefficients$mean)
  ms$pseudorsquared[i] <- round(allmodels[[i]]$pseudo.r.squared,4)
  ms$AICc[i] <- AICc(allmodels[[i]])
}

ms$deltaAICc <- ms$AICc - min(ms$AICc)
ms$wi <- exp(-0.5*ms$deltaAICc)/sum(exp(-0.5*ms$deltaAICc))
ms$er <- max(ms$wi)/ms$wi
ms <- ms[order(ms$modelnum),]
ms <- ms[order(ms$AICc),]   


relimp <- data.frame(variable = as.character(rep(NA,length(Cols))), wij = as.numeric(rep(NA,length(Cols))), stringsAsFactors = F)
for (i in 1:length(Cols)) {
  relimp$variable[i] <- Cols[i]
  relimp$wij[i] <- sum(ms[grep(Cols[i],ms$model),'wi'])
}
relimp <- relimp[rev(order(relimp$wij)),]

## plot model output ----
bestmodel <- allmodels[[ms$modelnum[1]]]
bestmodel2 <- allmodels[[ms$modelnum[2]]]
bestmodel3 <- allmodels[[ms$modelnum[3]]]
bestmodel4 <- allmodels[[ms$modelnum[4]]]
bestmodel5 <- allmodels[[ms$modelnum[5]]]

fecun$pred1 <- predict(bestmodel,newdata=fecun,type='response')
fecun$pred2 <- predict(bestmodel2,newdata=fecun,type='response')
fecun$pred3 <- predict(bestmodel3,newdata=fecun,type='response')
fecun$pred4 <- predict(bestmodel4,newdata=fecun,type='response')
fecun$pred5 <- predict(bestmodel5,newdata=fecun,type='response')

fecunice <- rbind(data.frame(cohort_year = 1955:1980, fecrate = rep(NA, length = length(1955:1980)), totpop = rep(NA, length = length(1955:1980)), abrate = rep(NA, length = length(1955:1980)), meancond = rep(NA, length = length(1955:1980))), fecunice)
fecunice$pred1 <- predict(bestmodel,newdata=fecunice,type='response')
fecunice$pred2 <- predict(bestmodel2,newdata=fecunice,type='response')
fecunice$pred3 <- predict(bestmodel3,newdata=fecunice,type='response')
fecunice$pred4 <- predict(bestmodel4,newdata=fecunice,type='response')
fecunice$pred5 <- predict(bestmodel5,newdata=fecunice,type='response')
#fecun$predpop <- predict(mfec1,newdata=fecun,type='response')

bmr2 <- round(summary(bestmodel)$pseudo.r.squared,3)
bmform <- bestmodel$formula
bm2r2 <- round(summary(bestmodel2)$pseudo.r.squared,3)
bm2form <- bestmodel2$formula
bm3r2 <- round(summary(bestmodel3)$pseudo.r.squared,3)
bm3form <- bestmodel3$formula
bm4r2 <- round(summary(bestmodel4)$pseudo.r.squared,3)
bm4form <- bestmodel4$formula
bm5r2 <- round(summary(bestmodel5)$pseudo.r.squared,3)
bm5form <- bestmodel5$formula

bmlab <-  bquote(R[italic(p)]^2 == .(format(bmr2, digits = 3)))
bm2lab <-  bquote(R[italic(p)]^2 == .(format(bm2r2, digits = 2)))
bm3lab <-  bquote(R[italic(p)]^2 == .(format(bm3r2, digits = 3)))
bm4lab <-  bquote(R[italic(p)]^2 == .(format(bm4r2, digits = 3)))
bm5lab <-  bquote(R[italic(p)]^2 == .(format(bm5r2, digits = 3)))


fec1lab <-   expression(paste("Pregnancy rate ~ abortion rate, R"[italic(p)]^2,"=",0.821))
fec2lab <-   expression(paste("Pregnancy rate ~ abortion rate + condition, R"[italic(p)]^2,"=",0.839))
fec3lab <-   expression(paste("Pregnancy rate ~ Population size + abortion rate, R"[italic(p)]^2,"=",0.824))
fec4lab <-   expression(paste("Pregnancy rate ~ Population size + abortion rate + condition, R"[italic(p)]^2,"=",0.841))
fec5lab <-   expression(paste("Pregnancy rate ~ condition, R"[italic(p)]^2,"=",0.517))

library(RColorBrewer)
mycolours <- brewer.pal(n = 5,name = "Set1")
mycolours

fecunice <- merge(fecun[,c('cohort_year', 'feclb', 'fecub')], fecunice)
fecunice[which(fecunice$cohort_year < 1981), 'feclb'] <- NA
fecunice[which(fecunice$cohort_year < 1981), 'fecub'] <- NA

fy <- ggplot(fecunice, aes(cohort_year, fecrate))
fy <- fy + geom_point()
fy <- fy + geom_linerange(aes(ymin = feclb, ymax = fecub))
fy <- fy + theme_set(theme_cowplot())
fy <- fy + labs(x = 'Cohort year', y = 'Pregnancy rate')
fy <- fy + geom_line(data = fecunice, aes(cohort_year, pred1), col = 'red', size = 1.5)
fy <- fy + geom_line(data = fecunice, aes(cohort_year, pred5), col = 'blue', size = 1.5)
fy <- fy + geom_line(data = fecunice, aes(cohort_year, pred4), col = 'darkgreen', size = 1.5)
fy <- fy + ylim(0, 1)
fy <- fy + xlim(1955, 2015)
fy


par(mfrow=c(1,1))
par(bty='n')
par(mar=c(4.5,4,1.5,1))  # set margins(bottom,left,top,right)
with(fecunice,plot(cohort_year,fecrate,pch=16,type='o',lty=3,xlab='Cohort year',ylab='fecunicedity rate',xlim=c(1950,2015),yaxp=c(0,1,2),ylim=c(0,1)))
axis(1,at=seq(1955,2015,by=5),labels=F,tck=-.01,col="#0000ff00",col.ticks='black')
axis(1,at=c(2010,2015),labels=F,tck=-.01,col="black",col.ticks="#0000ff00")
with(fecunice,lines(cohort_year,pred1,lwd=2,col=mycolours[1],type= 'o',pch=16))
with(fecunice,lines(cohort_year,pred2,lwd=2,col=mycolours[2],type= 'o',pch=16))   
with(fecunice,lines(cohort_year,pred3,lwd=2,col=mycolours[3],type= 'o',pch=16))   
with(fecunice,lines(cohort_year,pred4,lwd=2,col=mycolours[4],type= 'o',pch=16))   
with(fecunice,lines(cohort_year,pred5,lwd=2,col=mycolours[5],type= 'o',pch=16))   
with(fecunice,lines(cohort_year,fecrate,lty=3,type= 'o',pch=16))      
legend(1950, 0.4,cex=0.6, legend=c(fec1lab,fec2lab,fec3lab, fec4lab, fec5lab),bty='n',col=mycolours,lty=1,lwd=2)



par(mfrow=c(1,1))
par(bty='n')
par(mar=c(4.5,4,1.5,1))  # set margins(bottom,left,top,right)
with(fecunice,plot(cohort_year,fecrate,pch=16,type='o',lty=3,xlab='Cohort year',ylab='fecunicedity rate',xlim=c(1950,2015),yaxp=c(0,1,2),ylim=c(0,1)))
axis(1,at=seq(1955,2015,by=5),labels=F,tck=-.01,col="#0000ff00",col.ticks='black')
axis(1,at=c(2010,2015),labels=F,tck=-.01,col="black",col.ticks="#0000ff00")

with(fecunice,lines(cohort_year,pred5,lwd=2,col=mycolours[5],type= 'o',pch=16))   
with(fecunice,lines(cohort_year,fecrate,lty=3,type= 'o',pch=16))      
legend(1950, 0.4,cex=0.6, legend=c(fec1lab,fec2lab,fec3lab, fec4lab, fec5lab),bty='n',col=mycolours,lty=1,lwd=2)



