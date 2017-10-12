## load libraries ----
library(scales)
library(dplyr)
library(cowplot)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')
load('Rdata/data.Rdata')

data$doy <- yday(data$dates)

#data$doy <- ifelse(data$doy > 270, data$doy - 366, data$doy)


gamdata <- na.omit(data[,c('idsex', 'weight', 'ageclass', 'doy', 'sex','year','month', 'mdates', 'trueage')])
gamdata <- gamdata[which(!(gamdata$ageclass == '0' & gamdata$weight > 100)),]
gamdata <- gamdata[which(!(gamdata$ageclass == '1' & gamdata$weight > 120)),]
gamdata <- gamdata[which(!(gamdata$ageclass == '2' & gamdata$weight > 120)),]
gamdata <- gamdata[which(!(gamdata$ageclass == '3' & gamdata$weight > 150)),]

gamdata$stage <- ifelse(gamdata$trueage < 5, 'J','A')
gamdata$stage <- ifelse(gamdata$trueage == 0, 'YOY', gamdata$stage)

gamdata$decade <- ifelse(gamdata$year < 1990, '1980s',
                         ifelse(between(gamdata$year, 1990, 1999), '1990s',
                                ifelse(between(gamdata$year, 2000, 2009), '2000s', '2010s')))

gamdatayoy <- gamdata[which(gamdata$ageclass == '0'),]

juv <- 1:4
ad <- c(5:7, '8+')
gamdatajuv <- gamdata[which(gamdata$ageclass %in% juv),]
gamdataad <- gamdata[which(gamdata$ageclass %in% ad),]




# gamdata$age <- gamdata$ageclass
# gamdata[which(gamdata$age == '8+'), 'age'] <- 10
# gamdata$age <- as.numeric(gamdata$age)


library(mgcv)


mygamyoy <- gam(weight ~  s(doy,  bs = "cc" ), family = gaussian, data = gamdatayoy)
lmyoy <- lm(weight ~ doy, data = gamdatayoy)
mygamjuv <- gam(weight ~  s(doy,  bs = "cc" ) + as.factor(ageclass), family = gaussian, data = gamdatajuv)
mygamad <- gam(weight ~  s(doy,  bs = "cc" ) + as.factor(ageclass), family = gaussian, data = gamdataad)

# mygam1 <- gam(weight ~  s(doy,  bs = "cc" ) + s(doy, by =  trueage,  bs = "cc"), family = gaussian, data = gamdata)
# mygam12 <- gam(weight ~  te(doy,  trueage, bs = "cc" ), family = gaussian, data = gamdata)
# mygam2 <- gam(weight ~  s(doy,  bs = "cc" ) + as.factor(ageclass), family = gaussian, data = gamdata)
# 
# mygam3 <- gam(weight ~  s(doy,  bs = "cc" ) + as.factor(stage), family = gaussian, data = gamdata)


pdatyoy <- expand.grid(doy = min(gamdatayoy$doy):max(gamdatayoy$doy), ageclass = as.factor(unique(gamdatayoy$ageclass)))
pdatjuv <- expand.grid(doy = min(gamdata$doy):max(gamdata$doy), ageclass = as.factor(unique(gamdatajuv$ageclass)))
pdatad <- expand.grid(doy = min(gamdata$doy):max(gamdata$doy), ageclass = as.factor(unique(gamdataad$ageclass)))




pdatyoy$pred  <-  predict(mygamyoy,  newdata = pdatyoy, type = "response", se.fit = TRUE)[[1]]
pdatjuv$pred  <-  predict(mygamjuv,  newdata = pdatjuv, type = "response", se.fit = TRUE)[[1]]
pdatad$pred   <-  predict(mygamad,   newdata = pdatad , type = "response", se.fit = TRUE)[[1]]
#pdat$pred3  <-  predict(mygam3,  newdata = pdat, type = "response", se.fit = TRUE)[[1]]

pdat <- rbind(pdatyoy,pdatjuv, pdatad)
pdat <- (merge(pdat,unique(gamdata[c('doy', 'mdates')])))

datagroup <- data %>%
  select(ageclass, month, weight) %>%
  group_by(ageclass, month) %>%
  filter(complete.cases(ageclass)) %>%
  filter(month < 9000) %>%
  summarize(mw = mean(weight, na.rm = T))

datagroup$mdates <- (ymd(paste('2800',"-",datagroup$month,"-",15,sep = "")))

for (i in 1:nrow(datagroup)) {
  #if(!is.na(datagroup[i, 'dates'])) {
    if(datagroup$month[i] > 2) {
      datagroup$mdates[i] <- datagroup$mdates[i] - years(1)
   # }
  }
}



mw <- gamdata %>%
  group_by(decade, month, ageclass) %>%
  summarize(mw = mean(weight))

mw$mdates <- dmy(paste('15/', mw$month, '/2800'))
for (i in 1:nrow(mw)) {
  #if(!is.na(datagroup[i, 'dates'])) {
  if(mw$month[i] > 2) {
    mw$mdates[i] <- mw$mdates[i] - years(1)
    # }
  }
}


p <- ggplot(gamdata , aes(mdates, weight, colour = decade) ) 
p <- p + geom_point(alpha = 0.3)
p <- p + geom_line(data = mw, aes(mdates, mw, colour = decade) )
#p <- p + geom_line(data = pdat, aes(mdates, pred), size = 1.5, col = 'red')
#p <- p + geom_point(data = datagroup, aes(mdates, mw), size = 1.5, col = 'blue')
#p <- p + geom_line(data = pdat, aes(mdates, pred12, col = 'blue'))
#p <- p + theme_bw() + theme(legend.position="none")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_x_date(date_breaks = "1 months",labels = date_format("%b")) 
p <- p + xlab("Month")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p <- p + facet_wrap(~ageclass, ncol = 2)#, scales = "free_y")
print(p)

#save_plot("output/weight_gain_by_decade.pdf", p, base_height = 8, base_aspect_ratio = 1.3) # make room for figure legend)



mygamyoysd <- gam(weight ~  s(doy,  bs = "cc" ) + factor(decade) ,family=gaussian(link=log), data = gamdatayoy)
mygamjsd <- gam(weight ~  s(doy,  bs = "cc" ) + factor(decade) ,family=Gamma(link=log), data = gamdatajuv)
mygamasd <- gam(weight ~  s(doy,  bs = "cc" ) + factor(decade) , family=gaussian(link=log), data = gamdataad)

pdatyoy <- expand.grid(doy = min(gamdatayoy$doy):max(gamdatayoy$doy), stage = as.factor(unique(gamdatayoy$stage)), decade = as.factor(unique(gamdatayoy$decade)))
pdatjuv <- expand.grid(doy = min(gamdata$doy):max(gamdata$doy), stage = as.factor(unique(gamdatajuv$stage)), decade = as.factor(unique(gamdatayoy$decade)))
pdatad <- expand.grid(doy = min(gamdata$doy):max(gamdata$doy), stage = as.factor(unique(gamdataad$stage)), decade = as.factor(unique(gamdatayoy$decade)))

pdatyoy$pred  <-  predict(mygamyoysd,  newdata = pdatyoy, type = "response", se.fit = TRUE)[[1]]
pdatjuv$pred  <-  predict(mygamjsd,  newdata = pdatjuv, type = "response", se.fit = TRUE)[[1]]
pdatad$pred   <-  predict(mygamasd,   newdata = pdatad , type = "response", se.fit = TRUE)[[1]]

pdatsd <- rbind(pdatyoy,pdatjuv, pdatad)
pdatsd <- (merge(pdatsd,unique(gamdata[c('doy', 'mdates')])))


p <- ggplot(gamdata , aes(mdates, weight, colour = decade) ) 
p <- p + geom_point(alpha = 0.3)
p <- p + geom_line(data = pdatsd, aes(mdates, pred, colour = decade))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_x_date(date_breaks = "1 months",labels = date_format("%b")) 
p <- p + xlab("Month")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p <- p + facet_wrap(~stage, ncol = 1)#, scales = "free_y")
print(p)

save_plot("output/weight_gain_by_decade_gams.png", p, base_height = 8, base_width = 8) # make room for figure legend)
