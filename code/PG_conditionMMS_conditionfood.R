
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

capelin <- read.csv('data/capelin_acoustic.csv')
names(capelin) <- c('cohort_year', 'capbiomass')
fecun <- read.csv('data/fecundity_2014-feb20.csv', header = T)

load('Rdata/adf.Rdata')

## lag 
fecun$Capt1<-c(NA,fecun$Cap[1:(length(fecun$Cap)-1)])
fecun$Capt2<-c(NA,NA,fecun$Cap[1:(length(fecun$Cap)-2)])
fecun$ArcCodt1<-c(NA,fecun$ArcCod[1:(length(fecun$ArcCod)-1)])
fecun$ArcCodt2<-c(NA,NA,fecun$ArcCod[1:(length(fecun$ArcCod)-2)])
fecun$Sandt1<-c(NA,fecun$Sand[1:(length(fecun$Sand)-1)])
fecun$Sandt2<-c(NA,NA,fecun$Sand[1:(length(fecun$Sand)-2)])

## bootstrap condition data by cohort year ----
bootcond <- data.frame(cohort_year = sort(unique(adf$cohyear)), lcc = rep(NA, length(unique(adf$cohyear))), ucc = rep(NA, length(unique(adf$cohyear))))
for (i in sort(unique(adf$cohyear))) {
  dd <- quantile(base::sample(adf$relcond[which(adf$cohyear == i)], size = 100000, replace = T), c(0.025, 0.975))
  bootcond[which(bootcond$cohort_year == i), 'lcc'] <- dd[1]
  bootcond[which(bootcond$cohort_year == i), 'ucc'] <- dd[2]
}
rm(i)


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


## merge mean condition and CIs  ----
meancond <- merge(meancond, bootcond)
# for those years when meancond = NA, make lcc and ucc = NA
meancond[which(is.na(meancond$meancond)),c('lcc', 'ucc')] <- NA
rm(bootcond)
meancond <- merge(meancond, capelin, all = T)

## plot condition ----
cy <- ggplot(meancond, aes(cohort_year, meancond))
cy <- cy + geom_point()
cy <- cy + geom_linerange(aes(ymin = lcc, ymax = ucc))
cy <- cy + theme_set(theme_cowplot())
cy <- cy + labs(x = 'Cohort year', y = 'Relative condition')
cy




with(subset(meancond, cohort_year > 1988 & capbiomass < 800), plot((capbiomass), meancond, pch = 16, xlim = c(0,1000)))
with(subset(meancond, cohort_year > 1988 & capbiomass < 800), plot(log(capbiomass), meancond, pch = 16))
with(meancond, cor.test(log(capbiomass), meancond))    
with(meancond, cor.test((capbiomass), ilmeancond))     


meancond$ilmeancond <- gtools::inv.logit(meancond$meancond)


mm <- merge(meancond, fecun, by = 'cohort_year', all = T)
pairs(c(mm[,c('meancond','ice.1y.jan','Capt1','ArcCodt1','Sandt1', 'totpop', 'capbiomass')]),upper.panel=points, lower.panel=panelcor,diag.panel=panel.hist)


mm <- na.omit(mm[,c('cohort_year', 'meancond', 'Capt1', 'ArcCodt1', 'ice.1y.jan', 'lcc', 'ucc', 'capbiomass')])

m1 <- gam(meancond ~ te(Capt1)  + ice.1y.jan, data = mm)
m2 <- gam(meancond ~ te(ArcCodt1)  + ice.1y.jan, data = mm)
m3 <- gam(meancond ~ te(ArcCodt1) , data = mm)
m4 <- gam(meancond ~ te(ice.1y.jan) , data = mm)
m5 <- gam(meancond ~ te(ice.1y.jan) + ArcCodt1 , data = mm)
m6 <- gam(meancond ~ te(log(capbiomass)) , data = mm)

mm$Em1 <- predict(m1, newdata = mm)
mm$Em2 <- predict(m2, newdata = mm)
mm$Em3 <- predict(m3, newdata = mm)
mm$Em4 <- predict(m4, newdata = mm)
mm$Em5 <- predict(m5, newdata = mm)

se <- predict( m1 , se = TRUE)$se.fit
mm$Em1lb <- mm$Em1 - 1.96 * se
mm$Em1ub <- mm$Em1 + 1.96 * se
se <- predict( m2 , se = TRUE)$se.fit
mm$Em2lb <- mm$Em2 - 1.96 * se
mm$Em2ub <- mm$Em2 + 1.96 * se
se <- predict( m3 , se = TRUE)$se.fit
mm$Em3lb <- mm$Em3 - 1.96 * se
mm$Em3ub <- mm$Em3 + 1.96 * se
se <- predict( m4 , se = TRUE)$se.fit
mm$Em4lb <- mm$Em4 - 1.96 * se
mm$Em4ub <- mm$Em4 + 1.96 * se
se <- predict( m5 , se = TRUE)$se.fit
mm$Em5lb <- mm$Em5 - 1.96 * se
mm$Em5ub <- mm$Em5 + 1.96 * se


## plot condition ----
ay <- ggplot(mm, aes(cohort_year, meancond))
ay <- ay + geom_ribbon(data = mm, aes(ymin = Em5lb, ymax = Em5ub), alpha = 0.2, fill = 'red')
ay <- ay + geom_ribbon(data = mm, aes(ymin = Em2lb, ymax = Em2ub), alpha = 0.2, fill = 'blue')
ay <- ay + geom_point()
ay <- ay + geom_linerange(aes(ymin = lcc, ymax = ucc))
ay <- ay + geom_line(aes(cohort_year, Em5), colour = 'red')
ay <- ay + geom_line(aes(cohort_year, Em2), colour = 'blue')
ay <- ay + theme_set(theme_cowplot())
ay <- ay + labs(x = 'Year', y = 'Condition')
ay <- ay + annotate("text", x = 1997, y = 0.6, label = "condition ~ gam(smooth(ice.1y.jan) + ArcCodt1)", hjust = 0, size = 3)
ay <- ay + annotate("text", x = 1997, y = 0.65, label = "condition ~ gam(ice.1y.jan + smooth(ArcCodt1)", hjust = 0, size = 3)
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.6, yend = 0.6, colour = "red")
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.65, yend = 0.65, colour = "blue")
save_plot("output/condition-models.png", ay, base_aspect_ratio = 1, base_height = 4, base_width = 6) # make room for figure legend)

