## ADB August 17, 2017
## Purpose of script: analyze growth of harp seals
##                    - seasonal growth
##                    - decadal changes in growth, acconting for seasonality

## Data: 
##    Morph data from qsel_growth in mm_db_queries.accdb
##    Age data from qsel_age in mm_db_queries.accdb


## load libraries ----
library(ggplot2)
library(lubridate)

## read data ----
setwd('D:/Buren_files/GitHub/PG_Growth/')
age <- read.csv('data/PG_age.csv', header = T, as.is = T)
growth <- read.csv('data/PG_growth.csv', header = T, as.is = T)

## fix records ----
# One seal is coded as age 3, but everything points to it being much older, therefore remove age
age[which(age$idsex == '20082059F'),'cohortage'] <- NA
# one seal has length = 2300, turn it into NA
growth[which(growth$length == 2300),'length'] <- NA
# at least 3 seals have the wrong weight, turn them to NA
ids <-  c("19940928F", "20042978F", "20062138F")
growth[which(growth$idsex %in% ids), 'weight'] <- NA


## calculate volume ----
growth$volume <- with(growth, length*(girth^2))
## volume of a cylinder V = pi r^2 L
##                      Circumference (girth) = 2 pi r
##                                        r = girth / (2 pi)
growth$radius <- growth$girth / (2 * pi)
growth$volume <- with(growth, length * pi * (radius^2))

## merge data sets ----
data <- merge(growth, age[,c('idsex', 'cohortage')], by = 'idsex', all.x = T)
# there is a duplicate record
data <- data[!duplicated(data$idsex),]

## wrangle data set ----
# Create dates
data$dayf <- ifelse(is.na(data$day), 1, data$day)
data$dates <- ymd(paste(data$year,"-",data$month,"-",data$dayf,sep = ""))
data$mdates <- ymd(paste('2800',"-",data$month,"-",data$dayf,sep = ""))
for (i in 1:nrow(data)) {
  if(!is.na(data[i, 'dates'])) {
     if(data$month[i] > 2) {
       data$mdates[i] <- data$mdates[i] - years(1)
     }
  }
}

# Eliminate foetus, stillborns and starvlings, and YOY
# There is one seal of unknown sex, also eliminate it
# age = 99 does not include any individual in the previous categories. age=99 is unknown
eliminate <- c(91:93,99)
data <- data[which(!data$cohortage %in% eliminate),]
data <- subset(data, sex != 'U')
subset(data, cohortage == 99)

## there can not be seals with cohort age = 0 prior to march ----
data[which((data$mdates %in% ymd('2799-08-01'):ymd('2800-02-29')) & (data$cohortage == 0)), 'cohortage'] <-1 

## calculate true age ----
## between Sep 1 and Feb 29, substract 1 from cohort age
## between March 1 and Aug 30, set trueage = cohortage
data$trueage <- data$cohortage
#age0 <- data[which(data$cohortage == 0),]
#older <- data[which(data$cohortage != 0),]
t1 <- ymd('2799-09-01')
t2 <- ymd('2800-02-29')
data[data$mdates %in% t1:t2, 'trueage'] <- data[data$mdates %in% t1:t2, 'cohortage'] - 1
#data <- rbind(age0, older)



# create age classes
data$ageclass <- ifelse(data$trueage > 7, '8+', data$trueage)
# find outliers in weight by age class
g <- ggplot(data[which(!is.na(data$ageclass)),], aes(ageclass, weight))
g <- g + geom_violin() 
print(g)

unique(data[which(!is.na(data$ageclass)),'cohortage'])


## plots ----
p <- ggplot(data , aes(as.POSIXct(mdates), weight, colour = sex)) 
p <- p + geom_point()
p <- p + scale_x_datetime(date_breaks = "1 months") 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + xlab("Date")
p <- p + facet_wrap(~ageclass, ncol = 2, scales = "free_y")
print(p)



p2 <- ggplot(data , aes(as.POSIXct(dates), weight, colour = sex)) 
p2 <- p2 + geom_point()
p2 <- p2 + scale_x_datetime(date_breaks = "1 months") 
p2 <- p2 + theme_bw()
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- p2 + xlab("Date")
print(p2)

p3 <- ggplot(data , aes(as.POSIXct(dates), weight, colour = sex)) 
p3 <- p3 + geom_point()
p3 <- p3 + scale_x_datetime(date_breaks = "1 months") 
p3 <- p3 + theme_bw()
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3 <- p3 + xlab("Date")
p3 <- p3 + facet_wrap(~ageclass, ncol = 1, scales = "free_y")
print(p3)


## check ----
tocheck <- subset(data, trueage == 0)
check <- tocheck[which(tocheck$weight > 75),]   
tocheck <- subset(data, trueage == 1)
check <- rbind(check, tocheck[which(tocheck$weight > 100),])
tocheck <- subset(data, trueage == 2)
check <- rbind(check,tocheck[which(tocheck$weight > 105),] )  
tocheck <- subset(data, trueage == 3)
check <- rbind(check,tocheck[which(tocheck$weight > 150),])
check$issue <- '2heavy4age'

with(data,plot(length,weight))
check <- rbind(check, cbind(subset(data, length > 140 & weight < 15), issue = '2light4length'))
check <- rbind(check, cbind(subset(data, length > 110 & length < 120 & weight > 75), issue = '2heavy4length'))
check <- rbind(check, cbind(subset(data, length > 100 & length < 105 & weight > 57), issue = '2heavy4length'))
check <- rbind(check, cbind(subset(data, length > 128 & length < 140 & weight > 100), issue = '2heavy4length'))
check <- rbind(check, cbind(subset(data, length > 139 & length < 141 & weight > 117), issue = '2heavy4length'))
check <- rbind(check, cbind(subset(data, length > 149 & length < 151 & weight > 145), issue = '2heavy4length'))

rm(tocheck)
#check[order( check$idsex),]
check <- check[order( check$idsex),c('idsex', 'nafo', 'year', 'month', 'day', 'weight', 'length', 'pelage', 'cohortage', 'issue')]

#write.csv(check, 'output/check_weights.csv', row.names = F)


## estimate growth parameters ----
dat <- na.omit(data[,c('idsex', 'length', 'weight', 'blubberdepth', 'volume' ,'trueage', 'sex')])
dat <- subset(dat, trueage > 0)
with(dat, plot(log(volume), log(weight)))
## there is a clear outlier - log(volume) much smaller than 10 - remove it
dat[which(log(dat$volume) < 8), 'volume'] <- NA
#dat <- subset(dat, sex == 'M')

lmvw <- lm(weight ~ volume, data = dat)
a <- lmvw$coefficients[1]
b <- lmvw$coefficients[2]
with(dat, plot(volume, weight))
lines(c(0,5e+06), a + (c(0,5e+06))*b, col = 'red')

## calculate relative condition ----
dat$predweight <- a + (dat$volume)*b
dat$relcond <- with(dat, weight/predweight)

## merge back ----
data <- merge(data, dat[,c('idsex', 'relcond')], by = 'idsex', all.x = T)

## plot comdition ----
p <- ggplot(data , aes(as.POSIXct(mdates), relcond, colour = sex)) 
p <- p + geom_point()
p <- p + scale_x_datetime(date_breaks = "1 months") 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + xlab("Date")
p <- p + facet_wrap(~ageclass, ncol = 1, scales = "free_y")
print(p)



library(forecast)
myar <- Arima(data$weight, order = c(0,0,300), include.drift = TRUE, xreg = data$ageclass)
data$predw <- fitted(myar)

qqnorm(myar$residuals)
qqline(myar$residuals)
acf(data$weight)

p <- ggplot(data , aes(as.POSIXct(mdates), weight, colour = sex)) 
p <- p + geom_point(aes(as.POSIXct(data$mdates),data$predw, colour = data$sex))
p <- p + scale_x_datetime(date_breaks = "1 months") 
p <- p + geom_line(aes())
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + xlab("Date")
p <- p + facet_wrap(~ageclass, ncol = 1, scales = "free_y")
print(p)