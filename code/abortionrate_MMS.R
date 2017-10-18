inpath<-'D:/Buren_files/MEGA/papersAle/harps_fecundity/analyses2014/input'
outpath<-'D:/Buren_files/MEGA/papersAle/harps_fecundity/analyses2014/review2'
infile<-'fecundity_2014-feb20.csv'
montfile <- 'montecarlo-weightsforbetaregs.csv'
naofile <- 'Hurrell_NAO.txt'
bootfile <- 'bootCI_50000.csv'        ## 

linkfn <- 'loglog'

graphics.off()
op <- par(no.readonly = TRUE) # the whole list of settable par's
library(reshape)
library(betareg)
library(stringr)
library(MuMIn)
library(gplots)



setwd(inpath)
source('panelcor.r')
source(file="HighstatLibV7.R")
source("MCMCSupportHighstatV2.R")
montweight <- read.csv(montfile,header=T)
montweight <- montweight[,c('resample.level','weight')]
fecun<-read.csv(infile,header=T)
NAO <- read.table(naofile, header=F,sep=' ',strip.white=T,skip=1)
names(NAO) <- c('cohort_year','NAO')
bootci <- read.csv(bootfile,header=T)
names(bootci) <- c('cohort_year', 'feclb','fecub','ablb','abub')

fecun$fecrate<-fecun$preg/fecun$mature
fecun$totpreg<-with(fecun,EP+preg)
fecun$abrate<-fecun$EP/fecun$totpreg
#fecun$abrateold<-fecun$EP/fecun$totpreg
fecun$fecratet1<-c(NA,fecun$fecrate[1:(length(fecun$fecrate)-1)])
fecun$abratet1<-c(NA,fecun$abrate[1:(length(fecun$abrate)-1)])
fecun$totfood<-with(fecun,ArcCod+Cap+Sand)
fecun$totfoodt1<-c(NA,fecun$totfood[1:(length(fecun$totfood)-1)])
fecun$totfoodt2<-c(NA,NA,fecun$totfood[1:(length(fecun$totfood)-2)])

fecun$Capt1<-c(NA,fecun$Cap[1:(length(fecun$Cap)-1)])
fecun$Capt2<-c(NA,NA,fecun$Cap[1:(length(fecun$Cap)-2)])
fecun$ArcCodt1<-c(NA,fecun$ArcCod[1:(length(fecun$ArcCod)-1)])
fecun$ArcCodt2<-c(NA,NA,fecun$ArcCod[1:(length(fecun$ArcCod)-2)])
fecun$Sandt1<-c(NA,fecun$Sand[1:(length(fecun$Sand)-1)])
fecun$Sandt2<-c(NA,NA,fecun$Sand[1:(length(fecun$Sand)-2)])

fecun <- merge(fecun,montweight,by.x= "totpreg", by.y="resample.level",all.x=T)
fecun <- merge(fecun,NAO)
fecun <- merge(fecun,bootci,all.x=TRUE)

fecunice <- fecun
fecunice <- subset(fecunice, cohort_year<2014 & cohort_year>1995)
fecunice <- fecunice[order(fecunice$cohort_year),]
fecunice$abrate<-replace(fecunice$abrate,fecunice$abrate==0,1e-6)




pairs(c(fecunice[,c('ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','NAO')]),upper.panel=points, lower.panel=panelcor,diag.panel=panel.hist)
fecunicet1 <- fecunice[,c('abrate','ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','NAO','cohort_year','weight')]

########################################################
#Data exploration
#Outliers
MyVar <- c('abrate','ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','NAO')
MyXVar <- c('ice.1y.jan','Capt1','ArcCodt1','Sandt1','abratet1','NAO')
Mydotplot(fecunicet1[,MyVar])

#Collinearity
pairs(fecunicet1[,MyVar], lower.panel = panel.cor)

corvif(fecunicet1[ ,MyXVar])


VIF <-  corvif(fecunicet1[ ,MyXVar])
#setwd(outpath)
#write.csv(VIF,'VIF_abrate_review2.csv', row.names=T) 
#setwd(inpath)

########################################
### This code constructs and fits all submodels starting from a full model - main effects only ####
########################################
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
   
   
   
   
#############
## models with total food instead of the 3 separate species
######

#fecunicet1tf <- fecunice[,c('abrate','ice.1y.jan','totfoodt1','abratet1','NAO','cohort_year','weight')]
#
#########################################
#### This code constructs and fits all submodels starting from a full model - main effects only ####
#########################################
### obtain column names
#Cols <- names(fecunicet1tf)
### exclude response variable, weights and year
#Cols <- Cols[! Cols %in% c('cohort_year','weight','abrate')]
#n <- length(Cols)
### obtain all combinations of explanatory variables
#id <- unlist(
#        lapply(1:n,
#              function(i)combn(1:n,i,simplify=F)
#        )
#      ,recursive=F)
#
### apply these combinations into a formula, with abrate as response variable
#Formulas <- sapply(id,function(i)
#              paste("abrate~",paste(Cols[i],collapse="+"))
#            )
### number of models to be fitted            
#nmodels <-  length(Formulas)
### fit all models
#allmodelst1tf <- lapply(Formulas,function(i)
#    betareg(as.formula(i),data=fecunicet1tf,link=linkfn))
#
### obtain statistics for all models, including the expression of each model    
#mst1tf <- data.frame(modelnum=as.numeric(1:nmodels))
#for (i in 1:nrow(mst1tf)){
#   mst1tf$model[i] <- as.character(allmodelst1tf[[i]]$formula[3])
#   mst1tf$N[i] <-   allmodelst1tf[[i]]$nobs
#   mst1tf$K[i] <- length(allmodelst1tf[[i]]$coefficients$mean)
#   mst1tf$pseudorsquared[i] <- round(allmodelst1tf[[i]]$pseudo.r.squared,4)
#   mst1tf$AICc[i] <- AICc(allmodelst1tf[[i]])
#   mst1tf$LH[i] <- (allmodelst1tf[[i]]$loglik)
#   }
#
#
#
### combine the models with total food nd the models with the 3 prey by themselves
#mst1 <- rbind(mst1, mst1tf[grep('totfoodt1',mst1tf$model),])  
#
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


##############
#### assess "pretending" variables
################
# models with potential pretending variables
pretmod <- mst1[which(mst1$er<10 & mst1$er>1),c('modelnum')]
# variable number that is potentially pretending
pretvar <- 3

pretvariables <- data.frame(modelnum=pretmod)
i <- 0
for (mod in pretmod){
   i <- i+1
   pretvariables$model[i] <- as.character(allmodelst1[[mod]]$formula[3])
   pretvariables$diffLH[i] <- (allmodelst1[[mst1[which(mst1$deltaAICc==0),'modelnum']]]$loglik)-allmodelst1[[mod]]$loglik
   pretvariables$parestimate[i] <- allmodelst1[[mod]]$coefficients$mean[pretvar+1]
   pretvariables$lconfint[i] <- confint(allmodelst1[[mod]])[1+pretvar,1] 
   pretvariables$uconfint[i] <- confint(allmodelst1[[mod]])[1+pretvar,2]  
   }   
    
setwd(outpath)
#write.csv(mst1,'abrate-modelselection_lh.csv',row.names=F)
#write.csv(relimp,'abrate-relimp-variables.csv',row.names=F)
#write.csv(pretvariables,'pretending-variables.csv',row.names=F)

 

#
#
#fecunicet2 <- na.omit(fecunice[,c('abrate','ice.1y.jan','Capt2','ArcCodt2','Sandt2','abratet1','NAO','cohort_year','weight')])
#########################################
#### This code constructs and fits all submodels starting from a full model - main effects only ####
#########################################
### obtain column names
#Cols <- names(fecunicet2)
### exclude response variable, weights and year
#Cols <- Cols[! Cols %in% c('cohort_year','weight','abrate')]
#n <- length(Cols)
### obtain all combinations of explanatory variables
#id <- unlist(
#        lapply(1:n,
#              function(i)combn(1:n,i,simplify=F)
#        )
#      ,recursive=F)
#
### apply these combinations into a formula, with abrate as response variable
#Formulas <- sapply(id,function(i)
#              paste("abrate~",paste(Cols[i],collapse="+"))
#            )
### number of models to be fitted                        
#nmodels <-  length(Formulas)
### fit all models
#allmodelst2 <- lapply(Formulas,function(i)
#    betareg(as.formula(i),data=fecunicet2,link=linkfn))
#
### obtain statistics for all models, including the expression of each model        
#mst2 <- data.frame(modelnum=as.numeric(1:nmodels))
#for (i in 1:nrow(mst2)){
#   mst2$model[i] <- as.character(allmodelst2[[i]]$formula[3])
#   mst2$N[i] <-   allmodelst2[[i]]$nobs
#   mst2$K[i] <- length(allmodelst2[[i]]$coefficients$mean)
#   mst2$pseudorsquared[i] <- round(allmodelst2[[i]]$pseudo.r.squared,4)
#   mst2$AICc[i] <- AICc(allmodelst2[[i]])
#   }
#
#mst2$deltaAICc <- mst2$AICc - min(mst2$AICc)
#mst2$wi<-exp(-0.5*mst2$deltaAICc)/sum(exp(-0.5*mst2$deltaAICc))
#mst2$er <- max(mst2$wi)/mst2$wi
#mst2<-mst2[order(mst2$modelnum),]
#mst2<-mst2[order(mst2$AICc),]   
#
#

bestmodel <- allmodelst1[[mst1$modelnum[1]]]
bestmodel2 <- allmodelst1[[mst1$modelnum[2]]]
bestmodel3 <- allmodelst1[[mst1$modelnum[3]]]
bestmodel4 <- allmodelst1[[mst1$modelnum[4]]]
bestmodel5 <- allmodelst1[[mst1$modelnum[5]]]

fecun$pred1 <- predict(bestmodel,newdata=fecun,type='response')
fecun$pred2 <- predict(bestmodel2,newdata=fecun,type='response')
fecun$pred3 <- predict(bestmodel3,newdata=fecun,type='response')
fecun$pred4 <- predict(bestmodel4,newdata=fecun,type='response')
fecun$pred5 <- predict(bestmodel5,newdata=fecun,type='response')

#bmform <-   expression(paste("Rate of premature births ~ Ice coverage + Capelin",", R"[italic(p)]^2,"=",0.71,",  ",Delta*AIC[c]," = 0"))
#bm2form <-  expression(paste("Rate of premature births ~ Ice coverage + Capelin"," + Arctic cod"[(t-1)],", R"[italic(p)]^2,"=",0.77,",  ",Delta*AIC[c]," = 2.08"))
#bm3form <-  expression(paste("Rate of premature births ~ Ice coverage + Capelin"," + Sandlance"[(t-1)],", R"[italic(p)]^2,"=",0.7,",  ",Delta*AIC[c]," = 3.5"))
#bm4form <-  expression(paste("Rate of premature births ~ Ice coverage + Capelin"," + Abortion rate"[(t-1)],", R"[italic(p)]^2,"=",0.72,",  ",Delta*AIC[c]," = 3.62"))
#bm5form <-  expression(paste("Rate of premature births ~ Ice coverage + Capelin"," + NAO, R"[italic(p)]^2,"=",0.7,",  ",Delta*AIC[c]," = 3.87"))


# Special symbol plus text plus evaluated statements? Use bquote(). Put a ~ in between each component. Equal signs go in quotes. Wrap statements you want evaluated in .()

#main=bquote(mu[sepal ~ length] ~ "=" ~ .(mean(iris$Sepal.Length)))
#bmform  <- bquote(expression("Rate of premature births ~ Ice coverage + Capelin",", R"([italic(p)]^2),  "=") ~ .(summary(bestmodel)$pseudo.r.squared))

bmform  <- expression(paste("Abortion rate ~ Mid-winter ice + Capelin",", R"[italic(p)]^2,"=",0.71))
bm2form <- expression(paste("Abortion rate ~ Ice coverage + Capelin"," + Abortion rate"[(t-1)],", R"[italic(p)]^2,"=",0.69))
bm3form <- expression(paste("Abortion rate ~ Ice coverage + Capelin"," + Arctic cod"[(t-1)],", R"[italic(p)]^2,"=",0.69))
bm4form <- expression(paste("Abortion rate ~ Ice coverage + Capelin"," + Sandlance"[(t-1)],", R"[italic(p)]^2,"=",0.64)) 
bm5form <- expression(paste("Abortion rate ~ Ice coverage + Capelin"," + NAO, R"[italic(p)]^2,"=",0.64))

library(RColorBrewer)
mycolours<-brewer.pal(n=5,name="Set1")
mycolours


#B. Create X matrix 
tmpfec <- fecun[,c('abrate','ice.1y.jan','Capt1')]
tmpfec$abrate[is.na(tmpfec$abrate)] <- 0
tmpfec$ice.1y.jan[is.na(tmpfec$ice.1y.jan)] <- 0
tmpfec$Capt1[is.na(tmpfec$Capt1)] <- 0
Xmat <- model.matrix(~ ice.1y.jan+ Capt1, data = tmpfec)


#C. Calculate predicted values
eta <-  Xmat %*% coef(bestmodel)[1:length(coef(bestmodel))-1]

#D. Calculate standard errors (SE) for predicted values
#SE of fitted values are given by the square root of
#the diagonal elements of: X * cov(betas) * t(X)
SE <- sqrt(  diag(Xmat %*%vcov(bestmodel)[1:length(coef(bestmodel))-1,1:length(coef(bestmodel))-1] %*% t(Xmat))  )

#Glue them all together
#fecun$Pred <- fecun$Pred
etaplus <- eta + 1.96 * SE
etaminus <- eta - 1.96 * SE
fecun$pred1ub <-  exp(-exp(-etaplus))[,1]
fecun$pred1lb <-  exp(-exp(-etaminus))[,1]




abplotdata <-    subset(fecun,cohort_year>1995 & cohort_year<2014)

a <- ggplot (data = abplotdata, aes(x = cohort_year, abrate))
a <- a + geom_linerange(aes(ymin = ablb, ymax = abub))
a <- a + labs(x = 'Year', y = 'Abortion rate')
a <- a + geom_line(aes(cohort_year, pred1), colour =  'red', size = 1.5)
a <- a + geom_ribbon(data = abplotdata, aes(ymin = pred1lb, ymax = pred1ub), alpha = 0.2, fill = 'red')
a <- a + annotate("text", x = 1996.5 , y = 0.55, label = bmform, hjust = 0, size = 3)
a <- a + annotate("segment", x = 1995 , xend = 1996.3, y = 0.55, yend = 0.55, colour = 'red', size = 1.2)
a <- a + scale_x_continuous(breaks = seq(1995, 2015, 5), limits = c(1995, 2015))
a <- a + geom_point()
a
setwd('D:/Buren_files/GitHub/PG_Growth/')
save_plot("output/toGarry/abortion-as-paper_white.png", a, base_width = 8, base_height = 6, bg = 'white') # make room for figure legend)
save_plot("output/toGarry/abortion-as-paper_trans.png", a, base_width = 8, base_height = 6, bg = 'transparent') # make room for figure legend)



