
## load libraries ----
#library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(cowplot)
library(MuMIn)
library(betareg)
library(mgcv)
library(ggplot2)
library(RColorBrewer)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')
source("functions/HighstatLibV7.R")
source('functions/panelcor.r')
source("functions/MCMCSupportHighstatV2.R")

capelin <- read.csv('data/capelin_acoustic.csv')
names(capelin) <- c('cohort_year', 'capbiomass')
fecun <- read.csv('data/fecundity_2014-feb20.csv', header = T)

load('Rdata/adf.Rdata')
load('Rdata/diet.Rdata')

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
#meancond <- merge(meancond, capelin, all = T)

## filter diet data ----
prey <- c('arctic cod', 'capelin', 'mysid')
diet <- percbio %>%
  select(preycat, year, percbio) %>%
  filter(preycat %in% prey) %>%
  spread(preycat, percbio)
names(diet) <- c('cohort_year', 'darcticcod', 'dcapelin', 'dmysid')

## merge mean condition and diet  ----
meancond <- merge(meancond, diet, all = T)
rm(prey)

## merge with ice data ----
meancond <- merge(meancond, fecun[,c('cohort_year', 'ice.1y.jan')], all = T)

## plot condition ----
cy <- ggplot(meancond, aes(cohort_year, meancond))
cy <- cy + geom_point()
cy <- cy + geom_linerange(aes(ymin = lcc, ymax = ucc))
cy <- cy + theme_set(theme_cowplot())
cy <- cy + labs(x = 'Cohort year', y = 'Relative condition')
cy


mm <- meancond %>%
  filter(complete.cases(meancond))

pairs(c(mm[,c('meancond','darcticcod', 'dcapelin', 'dmysid', 'ice.1y.jan')]),upper.panel=points, lower.panel=panelcor,diag.panel=panel.hist)



m1 <- gam(meancond ~ te(dcapelin) , data = mm)
m2 <- gam(meancond ~ s(darcticcod) + s(ice.1y.jan, k = 5) , data = mm)
m3 <- gam(meancond ~ te(dmysid) , data = mm)
m4 <- gam(meancond ~ te(darcticcod)  + ice.1y.jan, data = mm)
m5 <- gam(meancond ~ s(darcticcod)  , data = mm)
m12 <- gam(meancond ~ s(ice.1y.jan), data = mm)






mm$Em1 <- predict(m1, newdata = mm)
mm$Em2 <- predict(m2, newdata = mm)
mm$Em3 <- predict(m3, newdata = mm)
mm$Em4 <- predict(m4, newdata = mm)
mm$Em5 <- predict(m5, newdata = mm)
mm$Em7 <- predict(m7, newdata = mm)
mm$Em8 <- predict(m8, newdata = mm)
mm$Em12 <- predict(m12, newdata = mm)

se <- predict( m1 , se = TRUE)$se.fit
mm$Em1lb <- mm$Em1 - 1.96 * se
mm$Em1ub <- mm$Em1 + 1.96 * se
se <- predict( m2 , se = TRUE)$se.fit
se <- c(0,se)
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
se <- predict( m7 , se = TRUE)$se.fit
mm$Em7lb <- mm$Em7 - 1.96 * se
mm$Em7ub <- mm$Em7 + 1.96 * se
se <- predict( m8 , se = TRUE)$se.fit
mm$Em8lb <- mm$Em8 - 1.96 * se
mm$Em8ub <- mm$Em8 + 1.96 * se
se <- predict( m12 , se = TRUE)$se.fit
mm$Em12lb <- mm$Em12 - 1.96 * se
mm$Em12ub <- mm$Em12 + 1.96 * se


mypalette <- c(brewer.pal(4,"Set1"))


## plot condition ----
ay <- ggplot(mm, aes(cohort_year, meancond))
ay <- ay + geom_text(aes(x = cohort_year, y = 1.5, label = paste('n = ',ncond, sep = ''), angle = 90), size = 3, hjust = 0.5)
# ay <- ay + geom_ribbon(data = mm, aes(ymin = Em5lb, ymax = Em5ub), alpha = 0.2, fill = mypalette[1])
ay <- ay + geom_ribbon(data = mm, aes(ymin = Em2lb, ymax = Em2ub), alpha = 0.2, fill = mypalette[1])
ay <- ay + geom_ribbon(data = mm, aes(ymin = Em7lb, ymax = Em7ub), alpha = 0.2, fill = mypalette[3])
ay <- ay + geom_ribbon(data = mm, aes(ymin = Em8lb, ymax = Em8ub), alpha = 0.2, fill =  mypalette[2])
ay <- ay + geom_point(aes(x = cohort_year, y = Em8), col =  mypalette[2], pch = 16, size = 0.6)
ay <- ay + geom_point(aes(x = cohort_year, y = Em7), col = mypalette[3], pch = 16, size = .6)
ay <- ay + geom_point(aes(x = cohort_year, y = Em12), col = mypalette[1], pch = 16, size = .6)
ay <- ay + geom_point(aes(x = cohort_year, y = Em2), col = mypalette[1], pch = 16, size = 1.3)
# ay <- ay + geom_line(aes(cohort_year, Em5), colour = mypalette[1])
ay <- ay + geom_line(aes(cohort_year, Em2), colour = mypalette[1])
ay <- ay + geom_line(aes(cohort_year, Em7), colour = mypalette[3])
ay <- ay + geom_line(aes(cohort_year, Em8), colour =  mypalette[2])
ay <- ay + geom_point()
ay <- ay + geom_linerange(aes(ymin = lcc, ymax = ucc))
ay <- ay + theme_set(theme_cowplot())
ay <- ay + labs(x = 'Year', y = 'Condition')
# ay <- ay + annotate("text", x = 1997, y = 0.6, label = "condition ~ gam(smooth(ice.1y.jan) + ArcCodt1)", hjust = 0, size = 3)
ay <- ay + annotate("text", x = 1997, y = 0.65, label = "condition ~ gam(smooth(ice.1y.jan)", hjust = 0, size = 2)
ay <- ay + annotate("text", x = 1997, y = 0.7, label = "condition ~ gam(smooth(ice.1y.jan) + smooth(ArcCodt1)", hjust = 0, size = 2)
ay <- ay + annotate("text", x = 1997, y = 0.75, label = "condition ~ gam(smooth(ice.1y.jan) + smooth(Capt1)", hjust = 0, size = 2)
# ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.6, yend = 0.6, colour = mypalette[1])
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.65, yend = 0.65, colour = mypalette[1])
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.7, yend = 0.7, colour = mypalette[3])
ay <- ay + annotate("segment", x = 1996, xend = 1996.8, y = 0.75, yend = 0.75, colour =  mypalette[2])
ay

save_plot("output/condition-models.png", ay, base_aspect_ratio = 1, base_height = 4, base_width = 6) # make room for figure legend)


## use data points instead of mean condition ----
cond <- adf[,c('cohyear', 'relcond')]
names(cond) <- c('cohort_year', 'relcond')
mm2 <- merge(cond, fecun, by = 'cohort_year', all = T)
mm2 <- mm2[order(mm2$cohort_year),]
mm2 <- (mm2[,c('cohort_year', 'relcond', 'Capt1', 'ArcCodt1', 'ice.1y.jan')])
mm2 <- merge(mm2, capelin, by = 'cohort_year', all = T)


mod1 <- gam(relcond ~ te(Capt1, ice.1y.jan, bs = 'tp'), data = mm2)
mod2 <- gam(relcond ~ s(Capt1)  + s(ice.1y.jan), data = mm2)
mod3 <- gam(relcond ~ te(ArcCodt1)  + s(ice.1y.jan), data = mm2)

xdat <- unique(mm2[,c('cohort_year', 'Capt1', 'ice.1y.jan')])



xdat$Em1 <- predict(mod1, newdata = xdat)
# mm$Em2 <- predict(m2, newdata = mm)
# mm$Em3 <- predict(m3, newdata = mm)
# mm$Em4 <- predict(m4, newdata = mm)
# mm$Em5 <- predict(m5, newdata = mm)
# mm$Em7 <- predict(m7, newdata = mm)
# mm$Em8 <- predict(m8, newdata = mm)



se <- predict( mod1 , se = TRUE)$se.fit
xdat$Em1lb <- xdat$Em1 - 1.96 * se
xdat$Em1ub <- xdat$Em1 + 1.96 * se



cp <- ggplot(mm2, aes(cohort_year, relcond))
#cp <- cp + geom_text(aes(x = cohort_year, y = 1.5, label = paste('n = ',ncond, sep = ''), angle = 90), size = 3, hjust = 0.5)
#cp <- cp + geom_ribbon(data = mm2, aes(ymin = Em1lb, ymax = Em1ub), alpha = 0.2, fill = mypalette[1])
# cp <- cp + geom_ribbon(data = mm, aes(ymin = Em2lb, ymax = Em2ub), alpha = 0.2, fill = mypalette[2])
# cp <- cp + geom_ribbon(data = mm, aes(ymin = Em7lb, ymax = Em7ub), alpha = 0.2, fill = mypalette[3])
# cp <- cp + geom_ribbon(data = mm, aes(ymin = Em8lb, ymax = Em8ub), alpha = 0.2, fill = mypalette[4])
cp <- cp + geom_point(data = xdat, aes(x = cohort_year, y = Em1), col = mypalette[1], pch = 16, size = 1.3)
# cp <- cp + geom_point(aes(x = cohort_year, y = Em7), col = mypalette[3], pch = 16, size = 1.3)
# cp <- cp + geom_point(aes(x = cohort_year, y = Em2), col = mypalette[2], pch = 16, size = 1.3)
# cp <- cp + geom_point(aes(x = cohort_year, y = Em5), col = mypalette[1], pch = 16, size = 1.3)
cp <- cp + geom_line(data = xdat, aes(cohort_year, Em1), colour = mypalette[1])
# cp <- cp + geom_line(aes(cohort_year, Em2), colour = mypalette[2])
# cp <- cp + geom_line(aes(cohort_year, Em7), colour = mypalette[3])
# cp <- cp + geom_line(aes(cohort_year, Em8), colour = mypalette[4])
cp <- cp + geom_point()
# cp <- cp + geom_linerange(aes(ymin = lcc, ymax = ucc))
cp <- cp + theme_set(theme_cowplot())
cp <- cp + labs(x = 'Year', y = 'Condition')
# cp <- cp + annotate("text", x = 1997, y = 0.6, label = "condition ~ gam(smooth(ice.1y.jan) + ArcCodt1)", hjust = 0, size = 3)
# cp <- cp + annotate("text", x = 1997, y = 0.65, label = "condition ~ gam(ice.1y.jan + smooth(ArcCodt1)", hjust = 0, size = 3)
# cp <- cp + annotate("text", x = 1997, y = 0.7, label = "condition ~ gam(smooth(ice.1y.jan) + smooth(ArcCodt1)", hjust = 0, size = 3)
# cp <- cp + annotate("text", x = 1997, y = 0.75, label = "condition ~ gam(smooth(ice.1y.jan) + smooth(Capt1)", hjust = 0, size = 3)
# cp <- cp + annotate("segment", x = 1996, xend = 1996.8, y = 0.6, yend = 0.6, colour = mypalette[1])
# cp <- cp + annotate("segment", x = 1996, xend = 1996.8, y = 0.65, yend = 0.65, colour = mypalette[2])
# cp <- cp + annotate("segment", x = 1996, xend = 1996.8, y = 0.7, yend = 0.7, colour = mypalette[3])
# cp <- cp + annotate("segment", x = 1996, xend = 1996.8, y = 0.75, yend = 0.75, colour = mypalette[4])
cp