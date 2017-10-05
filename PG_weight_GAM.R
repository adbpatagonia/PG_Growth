library(scales)
library(dplyr)

data$doy <- yday(data$dates)

#data$doy <- ifelse(data$doy > 270, data$doy - 366, data$doy)

gamdata <- na.omit(data[,c('idsex', 'weight', 'ageclass', 'doy', 'sex','year','month', 'mdates', 'trueage')])
gamdata <- gamdata[which(!(gamdata$ageclass == '0' & gamdata$weight > 100)),]
gamdata <- gamdata[which(!(gamdata$ageclass == '1' & gamdata$weight > 120)),]
gamdata <- gamdata[which(!(gamdata$ageclass == '2' & gamdata$weight > 120)),]
gamdata <- gamdata[which(!(gamdata$ageclass == '3' & gamdata$weight > 150)),]


gamdatayoy <- gamdata[which(gamdata$ageclass == '0'),]

juv <- 1:4
ad <- c(5:7, '8+')
gamdatajuv <- gamdata[which(gamdata$ageclass %in% juv),]
gamdataad <- gamdata[which(gamdata$ageclass %in% ad),]


gamdata$stage <- ifelse(gamdata$trueage < 5, 'J','A')

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


p <- ggplot(gamdata , aes(mdates, weight)) 
p <- p + geom_point()
p <- p + geom_line(data = pdat, aes(mdates, pred, col = 'red'), size = 1.5)
p <- p + geom_point(data = datagroup, aes(mdates, mw, col = 'blue'), size = 1.5)
#p <- p + geom_line(data = pdat, aes(mdates, pred12, col = 'blue'))
p <- p + theme_bw() + theme(legend.position="none")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_x_date(date_breaks = "1 months",labels = date_format("%b")) 
p <- p + xlab("Month")
p <- p + facet_wrap(~ageclass, ncol = 2)#, scales = "free_y")
print(p)


data2 <- merge(data, pdat[,c('doy','ageclass','pred')], by=c('doy','ageclass'), all.x = T)
data2 <- data2[!duplicated(data2$idsex),]

data2$resw <- data2$weight - data2$pred

