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

## merge data sets ----
data <- merge(growth, age[,c('idsex', 'cohortage')], by = 'idsex', all.x = T)

## wrangle data set ----
# Create dates
data$dayf <- ifelse(is.na(data$day), 1, data$day)
data$dates <- ymd(paste(data$year,"-",data$month,"-",data$dayf,sep = ""))
data$mdates <- ymd(paste('2800',"-",data$month,"-",data$dayf,sep = ""))
for (i in 1:nrow(data)) {
  if(!is.na(data[i, 'dates'])) {
     if(data$month[i] > 7) {
       data$mdates[i] <- data$mdates[i] - years(1)
     }
  }
}

# Eliminate foetus, stillborns and starvlings, and YOY
# There is one seal of unknown sex, also eliminate it
# age = 99 does not include any individual in the previous categories. age=99 is unknown
eliminate <- c(0, 91:93,99)
data <- data[which(!data$cohortage %in% eliminate),]
data <- subset(data, sex != 'U')
subset(data, cohortage == 99)

## calculate true age ----
## between Sep 1 and Feb 29, substract 1 from cohort age
## between March 1 and Aug 30, set trueage = cohortage
data$trueage <- data$cohortage
t1 <- ymd('2799-09-01')
t2 <- ymd('2800-02-29')
data[data$mdates %in% t1:t2, 'trueage'] <- data[data$mdates %in% t1:t2, 'cohortage'] - 1

# create age classes
data$ageclass <- ifelse(data$trueage > 7, '8+', data$trueage)
# find outliers in weight by age class
g <- ggplot(data[which(!is.na(data$ageclass)),], aes(ageclass, weight))
g <- g + geom_violin() 
print(g)

unique(data[which(!is.na(data$ageclass)),'cohortage'])


## There is an issue with cohort age! transform to real ages ----

## plots ----
p <- ggplot(data , aes(as.POSIXct(mdates), weight, colour = sex)) 
p <- p + geom_point()
p <- p + scale_x_datetime(date_breaks = "1 months") 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + xlab("Date")
p <- p + facet_wrap(~ageclass, ncol = 1, scales = "free_y")
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

write.csv(check, 'output/check_weights.csv', row.names = F)

