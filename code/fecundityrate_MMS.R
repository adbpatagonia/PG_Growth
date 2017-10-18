inpath<-'D:/Buren_files/MEGA/papersAle/harps_fecundity/analyses2014/input'
outpath<-'D:/Buren_files/MEGA/papersAle/harps_fecundity/analyses2014/review3'
infile<-'fecundity_2014-feb20.csv'
#montfile <- 'montecarlo-weightsforbetaregs.csv'
naofile <- 'Hurrell_NAO.txt'
bootfile <- 'bootCI_50000.csv'        ## 
#ppfile <- 'pup_production.csv'

linkfn <- 'log'

graphics.off()
op <- par(no.readonly = TRUE) # the whole list of settable par's
library(reshape)
library(betareg)
library(stringr)
library(MuMIn)
library(gplots)
library(ggplot2)


setwd(inpath)
source(file="HighstatLibV7.R")
source('panelcor.r')
source("MCMCSupportHighstatV2.R")
#montweight <- read.csv(montfile,header=T)
#montweight <- montweight[,c('resample.level','weight')]
fecun<-read.csv(infile,header=T)
NAO <- read.table(naofile, header=F,sep=' ',strip.white=T,skip=1)
names(NAO) <- c('cohort_year','NAO')
#pp <- read.csv(ppfile, header=T)
bootci <- read.csv(bootfile,header=T)
names(bootci) <- c('cohort_year', 'feclb','fecub','ablb','abub')


fecun$fecrate<-fecun$preg/fecun$mature
fecun$totpreg<-with(fecun,EP+preg)
fecun$abrate<-fecun$EP/fecun$totpreg
fecun$fecratet1<-c(NA,fecun$fecrate[1:(length(fecun$fecrate)-1)])
fecun$abratet1<-c(NA,fecun$abrate[1:(length(fecun$abrate)-1)])


#fecun <- merge(fecun,montweight,by.x= "mature", by.y="resample.level",all.x=T)
fecun <- merge(fecun,NAO)
fecun <- merge(fecun,bootci,all.x=TRUE)
#fecun$weight[which(is.na(fecun$weight))] <- max(montweight$weight)
fecun <- fecun[order(fecun$cohort_year),]
#fecun <- merge(fecun,pp,all.x=T)

fecunice <- fecun
fecunice <- subset(fecunice, cohort_year<2014 )

fecunice <- fecunice[,c('cohort_year','fecrate','totpop','abrate','NAO')]
fecunice <- na.omit(fecunice)
########################################################
#Data exploration
#Outliers
MyVar <- c('cohort_year','fecrate','totpop','abrate','NAO')
MyXVar <- c('totpop','abrate','NAO')
Mydotplot(fecunice[,MyVar])

#Collinearity
pairs(fecunice[,MyVar], lower.panel = panel.cor)

corvif(fecunice[ ,MyXVar])

VIF <-  corvif(fecunice[ ,MyXVar])
setwd('D:/Buren_files/MEGA/papersAle/harps_fecundity/analyses2014/review2')
#write.csv(VIF,'VIF_fecrate1.csv', row.names=T)
pairs(fecunice[ ,MyXVar], lower.panel = panel.cor)
setwd(inpath)
#
#
#### BEGIN MODEL FITS
#mfec1 <- betareg(fecrate~totpop,data=fecunice,link=linkfn)
#mfec2 <- betareg(fecrate~abrate+totpop,data=fecunice,link=linkfn)
#mfec3 <- betareg(fecrate~abrate,data=fecunice,link=linkfn)
#mfec4 <- betareg(fecrate~abrate+totpop+NAO,data=fecunice,link=linkfn)
#mfec5 <- betareg(fecrate~totpop+NAO,data=fecunice,link=linkfn)
#mfec6 <- betareg(fecrate~abrate+NAO,data=fecunice,link=linkfn)
#mfec7 <- betareg(fecrate~NAO,data=fecunice,link=linkfn)
#noweight.modfec1<-betareg(fecrate~totpop,data=fecunice,link=linkfn)
#noweight.modfec2<-betareg(fecrate~abrate+totpop,data=fecunice,link=linkfn)
#noweight.modfec9<-betareg(fecrate~abrate,data=fecunice,link=linkfn)
#
##mfec3<-betareg(fecrate~abrate+totpop+ice.1y.jan,data=fecunice,link=linkfn)
##mfec4<-betareg(fecrate~abrate+totpop+ice.1y.jan+fecratet1,data=fecunice,link=linkfn)
##mfec5<-betareg(fecrate~abrate+totpop+fecratet1,data=fecunice,link=linkfn)
##mfec6<-betareg(fecrate~totpop+ice.1y.jan,data=fecunice,link=linkfn)
##mfec7<-betareg(fecrate~totpop+ice.1y.jan+fecratet1,data=fecunice,link=linkfn)
##mfec8<-betareg(fecrate~totpop+fecratet1,data=fecunice,link=linkfn,control = betareg.control(start=c(-3.478808e-01,-2.458e-08,3.521e-01,41.28518 )))
##
##mfec10<-betareg(fecrate~ice.1y.jan,data=fecunice,link=linkfn)
##mfec11<-betareg(fecrate~fecratet1,data=fecunice,link=linkfn)
##mfec12<-betareg(fecrate~abrate+fecratet1,data=fecunice,link=linkfn)
##mfec13<-betareg(fecrate~ice.1y.jan+fecratet1,data=fecunice,link=linkfn,control = betareg.control(start=c(-2.99775392,-2.458e-08,3.521e-01,41.28518 )))
##mfec14<-betareg(fecrate~abrate+ice.1y.jan+fecratet1,data=fecunice,link=linkfn)
##mfec15<-betareg(fecrate~abrate+ice.1y.jan,data=fecunice,link=linkfn)
###mfec99<-betareg(fecrate~abrate+totpop,data=fecunice,link=linkfn, weights=mature)
##  
### END MODEL FITS
#


########################################
### This code constructs and fits all submodels starting from a full model - main effects only ####
########################################
## obtain column names
Cols <- names(fecunice)
## exclude response variable, weights and year
Cols <- Cols[! Cols %in% c('cohort_year','fecrate')]
n <- length(Cols)
## obtain all combinations of explanatory variables
id <- unlist(
        lapply(1:n,
              function(i)combn(1:n,i,simplify=F)
        )
      ,recursive=F)

## apply these combinations into a formula, with abrate as response variable
Formulas <- sapply(id,function(i)
              paste("fecrate~",paste(Cols[i],collapse="+"))
            )
## number of models to be fitted            
nmodels <-  length(Formulas)
## fit all models
allmodels <- lapply(Formulas,function(i)
    betareg(as.formula(i),data=fecunice,link=linkfn))

## obtain statistics for all models, including the expression of each model    
ms <- data.frame(modelnum=as.numeric(1:nmodels))
for (i in 1:nrow(ms)){
   ms$model[i] <- as.character(allmodels[[i]]$formula[3])
   ms$N[i] <-   allmodels[[i]]$nobs
   ms$K[i] <- length(allmodels[[i]]$coefficients$mean)
   ms$pseudorsquared[i] <- round(allmodels[[i]]$pseudo.r.squared,4)
   ms$AICc[i] <- AICc(allmodels[[i]])
   }

ms$deltaAICc <- ms$AICc - min(ms$AICc)
ms$wi<-exp(-0.5*ms$deltaAICc)/sum(exp(-0.5*ms$deltaAICc))
ms$er <- max(ms$wi)/ms$wi
ms<-ms[order(ms$modelnum),]
ms<-ms[order(ms$AICc),]   


 relimp <- data.frame(variable=as.character(rep(NA,length(Cols))), wij=as.numeric(rep(NA,length(Cols))), stringsAsFactors = F)
for (i in 1:length(Cols)){
   relimp$variable[i] <- Cols[i]
   relimp$wij[i] <- sum(ms[grep(Cols[i],ms$model),'wi'])
   }
relimp <- relimp[rev(order(relimp$wij)),]


#
#
#### CREATE TABLE FOR MODEL SELECTION
#library(MuMIn)
#ms <- data.frame(modelnum=as.numeric(str_sub(ls(pattern = "mfec"),5)))
#ms$model <- paste("mfec",ms$modelnum,sep="")
#for (i in 1:nrow(ms)){
#   ms$modelformula[i] <- as.character(eval(parse(text=ms$model[i]))$formula[3])
#   ms$N[i] <-   eval(parse(text=ms$model[i]))$nobs
#   ms$K[i] <- length(eval(parse(text=ms$model[i]))$coefficients$mean)
#   ms$pseudorsquared[i] <- round(summary(eval(parse(text=ms$model[i])))$pseudo.r.squared,4)
#   ms$LH[i] <- eval(parse(text=ms$model[i]))$loglik
#   ms$AICc[i] <- AICc(eval(parse(text=ms$model[i])))
#   }
#   
#ms$deltaAICc <- ms$AICc - min(ms$AICc)
#ms$wi<-exp(-0.5*ms$deltaAICc)/sum(exp(-0.5*ms$deltaAICc))
#ms$er <- max(ms$wi)/ms$wi
##
##
#ms<-ms[order(ms$modelnum),]
##
###write.csv(ms,'fecratemodels-step1-loglink-weight_var.csv',row.names=F)
#ms<-ms[order(ms$AICc),]
#### END CREATE TABLE FOR MODEL SELECTION
##
#### obtain column names
#Cols <- names(fecunice)
#### exclude response variable, weights and year
#Cols <- Cols[! Cols %in% c('cohort_year','weight','fecrate')]
#
#
# relimp <- data.frame(variable=as.character(rep(NA,length(Cols))), wij=as.numeric(rep(NA,length(Cols))), stringsAsFactors = F)
#for (i in 1:length(Cols)){
#   relimp$variable[i] <- Cols[i]
#   relimp$wij[i] <- sum(ms[grep(Cols[i],ms$modelformula),'wi'])
#   }
#relimp <- relimp[rev(order(relimp$wij)),]

 
#setwd(outpath)
#write.csv(relimp,'fecrate-step1-relimp-variables.csv',row.names=F)
#write.csv(ms,'fecrate-step1-modelselection.csv',row.names=F)

bestmodel <- allmodels[[ms$modelnum[1]]]
bestmodel3 <- allmodels[[ms$modelnum[3]]]

fecun$pred1 <- predict(bestmodel,newdata=fecun,type='response')
fecun$pred3 <- predict(bestmodel3,newdata=fecun,type='response')
fecunice$pred1 <- predict(bestmodel,newdata=fecunice,type='response')
fecunice$pred3 <- predict(bestmodel3,newdata=fecunice,type='response')
#fecun$predpop <- predict(mfec1,newdata=fecun,type='response')

bmr2 <- round(summary(bestmodel)$pseudo.r.squared,3)
bmform <- bestmodel$formula

bm3r2 <- round(summary(bestmodel3)$pseudo.r.squared,3)
bm3form <- bestmodel3$formula
bmlab <-  bquote(R[italic(p)]^2 == .(format(bmr2, digits = 3)))

bm3lab <-  bquote(R[italic(p)]^2 == .(format(bm3r2, digits = 3)))





fec1lab <-   expression(paste("Pregnancy rate ~ Population size + abortion rate, R"[italic(p)]^2,"=",0.79))

fec3lab <-   expression(paste("Pregnancy rate ~ Population size, R"[italic(p)]^2,"=",0.43))

library(RColorBrewer)
mycolours<-brewer.pal(n=3,name="Set1")
mycolours













#B. Create X matrix 
tmpfec <- fecun[,c('fecrate','abrate','totpop')]
tmpfec$abrate[is.na(tmpfec$abrate)] <- 0
tmpfec$totpop[is.na(tmpfec$totpop)] <- 0
Xmat <- model.matrix(~ totpop + abrate, data = tmpfec)


#C. Calculate predicted values
eta <-  Xmat %*% coef(bestmodel)[1:length(coef(bestmodel))-1]

#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
SE <- sqrt(  diag(Xmat %*%vcov(bestmodel)[1:length(coef(bestmodel))-1,1:length(coef(bestmodel))-1] %*% t(Xmat))  )

#Glue them all together
#fecun$Pred <- fecun$Pred
fecun$pred1ub <- exp(eta + 1.96 * SE)
fecun$pred1lb <- exp(eta - 1.96 * SE)
#####
#B. Create X matrix 
tmpfec <- fecun[,c('fecrate','abrate','totpop')]
tmpfec$abrate[is.na(tmpfec$abrate)] <- 0
tmpfec$totpop[is.na(tmpfec$totpop)] <- 0
Xmat <- model.matrix(~totpop, data = tmpfec)


#C. Calculate predicted values
eta <-  Xmat %*% coef(bestmodel3)[1:length(coef(bestmodel3))-1]

#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
SE <- sqrt(  diag(Xmat %*%vcov(bestmodel3)[1:length(coef(bestmodel3))-1,1:length(coef(bestmodel3))-1] %*% t(Xmat))  )

#Glue them all together
#fecun$Pred <- fecun$Pred
fecun$pred3ub <- exp(eta + 1.96 * SE)
fecun$pred3lb <- exp(eta - 1.96 * SE)



fecun$pred1lb <- ifelse(is.na(fecun$pred1),NA,fecun$pred1lb)
fecun$pred1ub <- ifelse(is.na(fecun$pred1),NA,fecun$pred1ub)
fecun$pred3lb <- ifelse(is.na(fecun$pred3),NA,fecun$pred3lb)
fecun$pred3ub <- ifelse(is.na(fecun$pred3),NA,fecun$pred3ub)


#setwd("D:/Buren_files/MEGA/papersAle/harps_fecundity/analyses2014/figs") 



setwd('D:/Buren_files/GitHub/PG_Growth/')

fplotdata <- subset(fecun, cohort_year < 2014)




fy <- ggplot(fplotdata, aes(cohort_year, fecrate))
fy <- fy + scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1))
fy <- fy + scale_x_continuous(breaks = seq(1950, 2015, 5), 
                              limits = c(1950, 2015),
                              labels = c(1950, rep("",1), 1960, rep("",1), 1970, rep("",1), 1980, rep("",1), 1990, rep("",1), 2000, rep("",1), 2010 , rep("",1))
                              )
fy <- fy + geom_point()
fy <- fy + geom_linerange(aes(ymin = feclb, ymax = fecub))
fy <- fy + theme_set(theme_cowplot())
fy <- fy + labs(x = 'Year', y = 'Pregnancy rate')
fy <- fy + geom_line(data = fplotdata, aes(cohort_year, pred1), col = 'red', size = 1.5)
fy <- fy + geom_line(data = fplotdata, aes(cohort_year, pred3), col = 'blue', size = 1.5, lty = 5 )
fy <- fy + geom_ribbon(data = fplotdata, aes(x = cohort_year, ymin = pred1lb, ymax = pred1ub), alpha = 0.2, fill = 'red')
fy <- fy + geom_ribbon(data = fplotdata, aes(x = cohort_year, ymin = pred3lb, ymax = pred3ub), alpha = 0.2, fill = 'blue')
fy <- fy + annotate("text", x = 1959.5 - 5  , y = 0.1, label = fec1lab, hjust = 0, size = 3)
fy <- fy + annotate("segment", x = 1955 - 5 , xend = 1959- 5, y = 0.1, yend = 0.1, colour = 'red', size = 1.2)
fy <- fy + annotate("text", x = 1959.5 - 5, y = 0.05, label = fec3lab, hjust = 0, size = 3)
fy <- fy + annotate("segment", x = 1955 - 5, xend = 1959 - 5, y = 0.05, yend = 0.05, colour = 'blue', size = 1.2, lty = 5)
fy


save_plot("output/toGarry/oldpreg-models_trans.png", fy, base_aspect_ratio = 1.4,  base_width = 6, bg = "transparent") # make room for figure legend)
save_plot("output/toGarry/oldpreg-models_white.png", fy, base_aspect_ratio = 1.4,  base_width = 6, bg = "white") # make room for figure legend)



