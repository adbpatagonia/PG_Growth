

## load libraries ----
#library(ggplot2)
library(lubridate)
library(dplyr)
library(cowplot)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')


age <- read.csv('data/PG_age.csv', header = T, as.is = T)
growth <- read.csv('data/PG_growth.csv', header = T, as.is = T)
maturity <-  read.csv('data/PG_maturity.csv', header = T, as.is = T)


## Mean monthly weight of foetuses (from database)
fwd <- 2.9 # december
fwj <- 6.3 # january
fwf <- 8.3 # february



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

library(cowplot)
w <- ggplot(adf, aes(exp(ll), exp(lw), colour = as.factor(month)))
w <- w + geom_point(alpha = 0.5)
w <- w + geom_line(data = pdata, aes(exp(ll), exp(lpw), colour = as.factor(month)), size = 1.5)
w <- w + theme_bw()
w <- w + labs(x = 'Length (cm)', y = 'Weight (kg)')
w <- w + scale_colour_discrete(name = "",
                             breaks = c("12", "1", "2"),
                             labels = c("December", "January", "February"))
w <- w + theme_set(theme_cowplot())
#save_plot("output/lw.png", w, base_aspect_ratio = 1.3, base_height = 6) # make room for figure legend)

lw <- ggplot(adf, aes((ll), (lw), colour = as.factor(month)))
lw <- lw + geom_point(alpha = 0.5)
lw <- lw + geom_line(data = pdata, aes((ll), (lpw), colour = as.factor(month)), size = 1.5)
lw <- lw + theme_bw()
lw <- lw + labs(x = 'log(Length (cm))', y = 'log(Weight (kg))')
lw <- lw + scale_colour_discrete(name = "",
                               breaks = c("12", "1", "2"),
                               labels = c("December", "January", "February"))
lw <- lw + theme_set(theme_cowplot())
#save_plot("output/llw.png", lw, base_aspect_ratio = 1.3, base_height = 6) # make room for figure legend)
save_plot("output/toGarry/llw-trans.png", lw, base_aspect_ratio = 1.3, base_height = 6, bg = 'transparent' ) # make room for figure legend)




## visualize condition ----
pregn <- 2:3
notpreg <- c(4:7,9)
ep <- 8
adf$preg <- ifelse(adf$fmaturity %in% pregn, 'A', 
                   ifelse(adf$fmaturity %in% notpreg, 'B',
                          ifelse(adf$fmaturity %in% ep, 'C', 'D')))

meanc <- adf %>%
  group_by(preg, cohyear) %>%
  summarize(meanc = mean(relcond))

labels <- c(A = "Pregnant", B = "Non pregnant", C = 'Early puppers')
cp <- ggplot(adf[which(adf$cohyear < 2015),], aes( factor(cohyear), relcond))
cp <- cp + geom_violin(trim = FALSE)
cp <- cp + geom_jitter(position = position_jitter(0.1)) + stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")
cp <- cp + xlab("Year") + ylab('Relative condition')
cp <- cp + geom_hline(yintercept = 1)
cp <- cp + facet_grid(preg ~ ., labeller=labeller(preg = labels))
cp <- cp + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
cp <- cp + theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust =0 ))
cp
save_plot("output/toGarry/condition_white.png", cp, base_width = 10, base_height = 6, bg = 'white') # make room for figure legend)



cp <- ggplot(adf[which(adf$cohyear < 2015),], aes( factor(cohyear), relcond))
cp <- cp + geom_violin(trim = FALSE)
cp <- cp + geom_jitter(position = position_jitter(0.1)) + stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")
cp <- cp + xlab("Year") + ylab('Relative condition')
cp <- cp + geom_hline(yintercept = 1)
cp <- cp + facet_grid(preg ~ ., labeller=labeller(preg = labels))
cp <- cp + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "transparent", colour = NA) )
cp <- cp + theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust =0 ))
cp
save_plot("output/toGarry/condition_trans.png", cp, base_width = 10, base_height = 6, bg = 'transparent') # make room for figure legend)



#save(adf, file = 'Rdata/adf.Rdata')
