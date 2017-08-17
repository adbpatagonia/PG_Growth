## ADB August 17, 2017
## Purpose of script: analyze growth of harp seals
##                    - seasonal growth
##                    - decadal changes in growth, acconting for seasonality

## Data: 
##    Morph data from qsel_growth in mm_db_queries.accdb
##    Age data from qsel_age in mm_db_queries.accdb


## load libraries
library(ggplot2)
library(lubridate)

## read data
setwd('D:/Buren_files/GitHub/PG_Growth/')
age <- read.csv('data/PG_age.csv', header = T, as.is=T)
growth <- read.csv('data/PG_growth.csv', header = T, as.is=T)

## merge data sets
data <- merge(growth, age[,c('idsex', 'cohortage')], by='idsex', all.x=T)

## wrangle data set
# Create dates
data$dayf <- ifelse(is.na(data$day), 1, data$day)
data$dates <- ymd(paste(data$year,"-",data$month,"-",data$dayf,sep=""))
data$mdates <- ymd(paste('3000',"-",data$month,"-",data$dayf,sep=""))

# mdate for 1 record failed to parse - assign it the value that it parsed for date
data$mdates[which(data$idsex == '19920095F')] <- data$dates[which(data$idsex == '19920095F')]

p <- ggplot(data , aes(as.POSIXct(mdates), weight)) 
p <- p + geom_point()
p <- p + scale_x_datetime(date_breaks="1 months") 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + xlab("Date")
p <- p + 
