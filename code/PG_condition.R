library(stringr)

## subset data ----
m <- c(1,2,12)
adf <- na.omit(data[which(data$cohortage >4 & data$sex =='F' & data$month %in% m & data$length > 130), c('idsex', 'month', 'length', 'weight', 'cohortage')])

#adf <- data %>%
 # filter(cohortage > 4, sex == 'F' & month %in% m)
#%>%
#  select( month, length, weight, cohortage) %>% 
  #filter(complete.cases(.))

adf$ll <- log(adf$length); adf$lw <- log(adf$weight)

adfd <- adf %>%
  filter(month == 12)

adfj <- adf %>%
  filter(month == 1)

adff <- adf %>%
  filter(month == 2)


## fit models and get predicted weight ----
fitd <- lm( lw ~ ll, data = adfd)
adfd$predw <- exp(predict(fitd, type = 'response'))
fitj <- lm( lw ~ ll, data = adfj)
adfj$predw <- exp(predict(fitj, type = 'response'))
fitf <- lm( lw ~ ll, data = adff)
adff$predw <- exp(predict(fitf, type = 'response'))

fitdp <- lm( lw ~ ll + I(ll^2)  , data = adfd)
adfd$predwp <- exp(predict(fitdp, type = 'response'))
fitjp <- lm( lw ~ ll + I(ll^2) , data = adfj)
adfj$predwp <- exp(predict(fitjp, type = 'response'))
fitfp <- lm( lw ~ ll + I(ll^2) , data = adff)
adff$predwp <- exp(predict(fitfp, type = 'response'))

fitds <- segmented(fitd, seg.Z = ~ll, psi = 5)
adfd$predws <- exp(predict(fitds, type = 'response'))
fitjs <- segmented(fitj, seg.Z = ~ll, psi = 5)
adfj$predws <- exp(predict(fitjs, type = 'response'))
fitfs <- segmented(fitf, seg.Z = ~ll, psi = 5)
adff$predws <- exp(predict(fitfs, type = 'response'))



fitdpol <- gamlss(lw ~ poly(ll,2), data = adfd, family=NO)
fitjpol <- gamlss(lw ~ poly(ll,2), data = adfj, family=NO)
fitfpol <- gamlss(lw ~ poly(ll,2), data = adff, family=NO)

adfd$predwp <- exp(predict(fitdpol, type = 'response'))
adfj$predwp <- exp(predict(fitjpol, type = 'response'))
adff$predwp <- exp(predict(fitfpol, type = 'response'))

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


pdatad$lpws <- predict(fitds, newdata = pdatad, type = 'response')  
pdataj$lpws <- predict(fitjs, newdata = pdataj, type = 'response')  
pdataf$lpws <- predict(fitfs, newdata = pdataf, type = 'response')  


pdatad$lpwp <- predict(fitdp, newdata = pdatad, type = 'response')  
pdataj$lpwp <- predict(fitjp, newdata = pdataj, type = 'response')  
pdataf$lpwp <- predict(fitfp, newdata = pdataf, type = 'response') 



pdata <- rbind(pdataj, pdataf, pdatad)
# pdatap <- rbind(pdatajp, pdatafp, pdatadp)
# pdatas <- rbind(pdatajs, pdatafs, pdatads)



w <- ggplot(adf, aes(exp(ll), exp(lw), colour = as.factor(month)))
w <- w + geom_point()
w <- w + geom_line(data= pdata, aes(exp(ll), exp(lpw), colour = as.factor(month)), size = 1.5, lty =2)
#w <- w + geom_line(data= pdata, aes(exp(ll), exp(lpws), colour = as.factor(month)), size = 1.5, lty =2)
w <- w + geom_line(data= pdata, aes(exp(ll), exp(lpwp), colour = as.factor(month)), size = 1.5, lty =1)
w <- w + theme_bw()
w

## visualize condition ----

cp <- ggplot(adf, aes( factor(cohyear), relcond))
cp <- cp + geom_violin()
cp <- cp + stat_summary(aes(factor(cohyear)), fun.y = mean, geom = "point", fill = "black", shape = 21, size = 1.5)
cp <- cp + xlab("Year") + ylab('Condition')
#cp <- cp + facet_grid(cohortage ~ ., scales = "free_y")
cp <- cp + theme_bw()
cp <- cp + geom_hline(yintercept = 1)
print(cp)

## CONCLUSION FROM EXPLORATIONS ----
# 1. STICK TO W = aL^b
# use mature females