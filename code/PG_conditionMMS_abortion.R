
## load libraries ----
#library(ggplot2)
library(lubridate)
library(dplyr)
library(cowplot)
library(MuMIn)
library(betareg)
library(mgcv)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')
source("functions/HighstatLibV7.R")
source('functions/panelcor.r')
source("functions/MCMCSupportHighstatV2.R")

load('Rdata/adf.Rdata')


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
ncond <-adf %>%
  group_by(cohyear) %>%
  group_size()
meancond <- cbind(meancond, ncond)
rm(ncond)
meancond[which(meancond$ncond < 10), 'meancond'] <- NA


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


rm(bootcond, meancond, bootci, adf)

## plot fecundity ----
ay <- ggplot(fecun, aes(cohort_year, abrate))
ay <- ay + geom_point()
ay <- ay + geom_linerange(aes(ymin = ablb, ymax = abub))
ay <- ay + theme_set(theme_cowplot())
ay <- ay + labs(x = 'Cohort year', y = 'Abortion rate')
ay


cy <- ggplot(fecun, aes(cohort_year, meancond))
cy <- cy + geom_point()
cy <- cy + geom_linerange(aes(ymin = lcc, ymax = ucc))
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



## model abortion as a betareg of condition ----
linkfn <- 'logit'
mfecun <- na.omit(fecun[,c('meancond', 'abrate')])
mfecun$mabrate <- mfecun$abrate + .000000001
betab <- betareg(mabrate ~ meancond, data = mfecun,link = linkfn)
xdat <- data.frame(meancond = seq(min(na.omit(fecun$meancond)), max(na.omit(fecun$meancond)), length.out =100))
xdat$Eabrate <- predict(betab, newdata = xdat, type = 'response')

ca <- ggplot(xdat, aes(meancond, Eabrate))
ca <- ca + geom_line(colour = 'red', size = 1.5)
ca <- ca + ylim(0, 0.6)
ca <- ca + xlim(0.8, 1.2)
ca <- ca + geom_point(data = fecun, aes(meancond, abrate))
ca <- ca + geom_linerange(data = fecun, aes(x = meancond, y = abrate, ymin = ablb, ymax = abub))
ca <- ca + labs(x = 'Relative condition', y = 'Abortion rate')
ca <- ca + theme_set(theme_cowplot())
ca
save_plot("output/abortion-condition-fit.png", ca, base_aspect_ratio = 1, base_height = 6) # make room for figure legend)

ca <- ggplot(fecun, aes(meancond, abrate))
ca <- ca + ylim(0, 0.6)
ca <- ca + xlim(0.8, 1.2)
ca <- ca + geom_point()
ca <- ca + geom_linerange(data = fecun, aes(x = meancond, y = abrate, ymin = ablb, ymax = abub))
ca <- ca + labs(x = 'Relative condition', y = 'Abortion rate')
ca <- ca + theme_set(theme_cowplot())
ca
save_plot("output/abortion-condition.png", ca, base_aspect_ratio = 1, base_height = 6) # make room for figure legend)



with(mfecun, plot(meancond, abrate,pch=16,ylim=c(0,.6),xlab='Mean relative condition',ylab='Abortion rate'))
#with(mfecun,plotCI(x=meancond,y=abrate,ui=abub,li=ablb,type='p',lty=1,sfrac=0,pch=16,gap=0,cex=0.1,add=T)) 
#   with(mfecun,lines(relcond,Eabs,col='red',lwd=2))
# with(mfecun,lines(meancond[I],Ebetabs[I],col='red',lwd=2))  
# with(mfecun,lines(meancond[I],Eabs[I],col='blue',lwd=2))  







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
fecunice$labrate <- gtools::logit(fecunice$abrate + 0.0000001)


pairs(c(fecunice[,c('ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1', 'meancond')]), #pch = 16, 
      upper.panel = points, lower.panel = panelcor,diag.panel=panel.hist)
fecunicet1 <- na.omit(fecunice[,c('abrate','ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','cohort_year', 'meancond')])

#Data exploration
#Outliers
MyVar <- c('abrate','ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','meancond')
MyXVar <- c('ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','meancond')
#Mydotplot(fecunicet1[,MyVar])

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


bestmodel <- allmodelst1[[mst1$modelnum[1]]]

plot(fecunicet1$abrate)
#points(fitted(bm), pch = 16, col = 'red')
points(fitted(bestmodel), pch = 16, col = 'purple')

## gams ----


## define link function
linkfn <- 'logit'
## obtain column names
Cols <- names(fecunicet1)
## exclude response variable, weights and year
Cols <- Cols[!Cols %in% c('cohort_year','weight','abrate', 'meancond')]
Cols <- c(Cols, 's(meancond)')
n <- length(Cols)
## obtain all combinations of explanatory variables
id <- unlist(
  lapply(1:n,
         function(i)combn(1:n,i,simplify = F)
  )
  ,recursive = F)

## apply these combinations into a formula, with abrate as response variable
Formulas <- sapply(id, function(i)
  paste("abrate ~ ", paste(Cols[i], collapse = "+"))
)
## number of models to be fitted            
nmodels <-  length(Formulas)
## fit all models
allgammodelst1 <- lapply(Formulas,function(i)
  gam(as.formula(i), data = fecunicet1))



## obtain statistics for all gam models, including the expression of each model    
mst1gam <- data.frame(modelnum = as.numeric(1:nmodels))
for (i in 1:nrow(mst1gam)) {
  mst1gam$model[i] <- as.character(allgammodelst1[[i]]$formula[3])
  mst1gam$N[i] <-   dim(allgammodelst1[[i]]$model)[1]
#  mst1gam$K[i] <- length(allgammodelst1[[i]]$coefficients$mean)
#  mst1gam$pseudorsquared[i] <- round(allgammodelst1[[i]]$pseudo.r.squared,4)
  mst1gam$AICc[i] <- AICc(allgammodelst1[[i]])
  mst1gam$LH[i] <- logLik(allgammodelst1[[i]])[1]
}

### compute AICc and other stats for all models 
mst1gam$deltaAICc <- mst1gam$AICc - min(mst1gam$AICc)

mst1gam$wi <- exp(-0.5*mst1gam$deltaAICc)/sum(exp(-0.5*mst1gam$deltaAICc))
mst1gam$er <- max(mst1gam$wi)/mst1gam$wi
mst1gam <- mst1gam[order(mst1gam$modelnum),]
mst1gam <- mst1gam[order(mst1gam$AICc),]  

mst1gam[which(mst1gam$er < 10), c('model','deltaAICc','er', 'LH')]


Cols <- c( names(fecunicet1))
Cols <- unique(Cols[! Cols %in% c('cohort_year','weight','abrate')])
relimpgam <- data.frame(variable=as.character(rep(NA,length(Cols))), wij=as.numeric(rep(NA,length(Cols))), stringsAsFactors = F)
for (i in 1:length(Cols)){
  relimpgam$variable[i] <- Cols[i]
  relimpgam$wij[i] <- sum(mst1gam[grep(Cols[i],mst1gam$model),'wi'])
}
relimpgam <- relimpgam[rev(order(relimpgam$wij)),]

bestmodel <- allgammodelst1[[mst1gam$modelnum[1]]]
bestmodel2 <- allmodelst1[[mst1$modelnum[3]]]
bestmodel2 <- betareg(abrate ~ Capt1 + ice.1y.jan, data = fecunicet1, link = 'loglog' )


fecun$pred1 <- predict(bestmodel,newdata=fecun,type='response')
fecun$pred2 <- predict(bestmodel2,newdata=fecun,type='response')



## obtain CIs for predictions ----
plotfecun <- subset(fecun,cohort_year > 1995 & cohort_year < 2012)
## for ice + capelin model
# Create X matrix 
tmpfec <- plotfecun[,c('abrate','Capt1','ice.1y.jan')]
Xmat <- model.matrix(~ Capt1 + ice.1y.jan, data = tmpfec)
# Calculate predicted values
eta <-  Xmat %*% coef(bestmodel2)[1:length(coef(bestmodel2))-1]
#Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
SE <- sqrt(  diag(Xmat %*%vcov(bestmodel2)[1:length(coef(bestmodel2))-1,1:length(coef(bestmodel2))-1] %*% t(Xmat))  )
#Glue them all together
etaplus <- eta + 1.96 * SE
etaminus <- eta - 1.96 * SE
plotfecun$pred2ub <-  exp(-exp(-etaplus))
plotfecun$pred2lb <-  exp(-exp(-etaminus))


# ## for s(meancond) model
se <- predict( bestmodel , se = TRUE)$se.fit
plotfecun$pred1lb <- plotfecun$pred1 - 1.96 * se
plotfecun$pred1ub <- plotfecun$pred1 + 1.96 * se


plotfecun[which(plotfecun$pred1lb < 0), 'pred1lb'] <- 0
plotfecun[which(plotfecun$pred2lb < 0), 'pred2lb'] <- 0






ay <- ggplot(plotfecun, aes(cohort_year, abrate))
ay <- ay + geom_point()
ay <- ay + geom_linerange(aes(ymin = ablb, ymax = abub))
ay <- ay + theme_set(theme_cowplot())
ay <- ay + labs(x = 'Year', y = 'Abortion rate')
ay <- ay + geom_line(data = plotfecun, aes(x = cohort_year, y = pred1), col = 'red')
ay <- ay + geom_line(data = plotfecun, aes(x = cohort_year, y = pred2), col = 'blue')
ay <- ay + geom_point(data = plotfecun, aes(x = cohort_year, y = pred2), col = 'blue', pch = 16, size = 1.3)
ay <- ay + geom_point(data = plotfecun, aes(x = cohort_year, y = pred1), col = 'red', pch = 16, size = 1.3)
ay <- ay + geom_ribbon(data = plotfecun, aes(ymin = pred2lb, ymax = pred2ub), alpha = 0.2, fill = 'blue')
ay <- ay + geom_ribbon(data = plotfecun, aes(ymin = pred1lb, ymax = pred1ub), alpha = 0.2, fill = 'red')
ay <- ay + annotate("text", x = 1997, y = 0.6, label = "Abortion rate ~ gam(condition)", hjust = 0, size = 2.6)
ay <- ay + annotate("text", x = 1997, y = 0.58, label = "Abortion rate ~ betareg(Capelin + January ice)", hjust = 0, size = 2.6)
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.6, yend = 0.6, colour = "red")
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.58, yend = 0.58, colour = "blue")

save_plot("output/abortion-models.png", ay, base_aspect_ratio = 1, base_height = 4, base_width = 6) # make room for figure legend)
