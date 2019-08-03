# Packages List
library(readr)
library(fUnitRoots)
library(tseries)
library(lmtest)
library(TSA)
library(forecast)
library(CombMSC)
library(fGarch)
library(rugarch)

# Functions List

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH")[1]){
  # If you have an output from arima() function use class = "ARIMA"
  # If you have an output from garch() function use class = "GARCH"
  # If you have an output from ugarchfit() function use class = "ARMA-GARCH"
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  acf(res.model,main="ACF of standardised residuals")
  pacf(res.model,main="PACF of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
}

MASE = function(observed , fitted ){
  # observed: Observed series on the forecast period
  # fitted: Forecast values by your model
  Y.t = observed
  n = length(fitted)
  e.t = Y.t - fitted
  sum = 0
  for (i in 2:n){
    sum = sum + abs(Y.t[i] - Y.t[i-1] )
  }
  q.t = e.t / (sum/(n-1))
  MASE = data.frame( MASE = mean(abs(q.t)))
  return(list(MASE = MASE))
}

# Reading the dataset
Bitcoin<- read_csv("C:/Users/Syed Hassan Afsar/Downloads/RMIT/2nd Semester/Time Series Analysis/Final Project/Bitcoin_Historical_Price.csv")
index <- seq(as.Date("2013-04-27"), as.Date("2019-02-24"), by="day")

# Convertiong the dataset to timeseries object
Bitcoin_ts <- ts(Bitcoin$Close, start = c(2013, as.numeric(format(index[1], "%j"))), frequency = 365)

class(Bitcoin_ts)

# Checking for missing values and imputing

which(is.na(Bitcoin_ts))

Bitcoin_ts[1773] <- 11573.30

Bitcoin_ts[1774] <- 10779.90

Bitcoin_ts[1775] <- 9965.57

Bitcoin_ts[1776] <- 9395.01

Bitcoin_ts[1777] <- 9337.55

Bitcoin_ts[1778] <- 8866.00

Bitcoin_ts[1779] <- 9578.63

Bitcoin_ts[1780] <- 9205.12

Bitcoin_ts[1781] <- 9194.85

Bitcoin_ts[1782] <- 8269.81

which(is.na(Bitcoin_ts))

# Visualisation using ACF, PACF plot

qqnorm(Bitcoin_ts,main="Q-Q Normal Plot of Bitcoin series")
qqline(Bitcoin_ts) # Fat tails is in accordance with volatiliy clustering

par(mfrow=c(1,2))
acf(Bitcoin_ts, main="The sample ACF plot for Bitcoin series")
pacf(Bitcoin_ts, main="The sample PACF plot for Bitcoin series")
eacf(Bitcoin_ts)


# Checking for ARCH effect

McLeod.Li.test(y=Bitcoin_ts,main="McLeod-Li Test Statistics for Bitcoin series")


# Transformations

BC.Bitcoin <- BoxCox.ar(Bitcoin_ts, method="yule-walker")

# 0 

# log transformation

par(mfrow=c(1,2))
acf(log(Bitcoin_ts))
pacf(log(Bitcoin_ts))

# stil there is a trend / decaying pattern
# high first lag in the pacf

qqnorm(log(Bitcoin_ts), main = "QQplot for natural log of bitcoin series")
qqline(log(Bitcoin_ts))

adf.test(log(Bitcoin_ts))
# 0.7348
# Still not stationary , hence we go for differencing

diff.log.bitcoin = diff(log(Bitcoin_ts),difference = 1)
plot(diff.log.bitcoin, type = 'o', ylab = 'Closing price')

order = ar(diff(diff.log.bitcoin))$order
adfTest(diff.log.bitcoin,lags = order, title = NULL, description = NULL) #Stationary with lag order 32

McLeod.Li.test(y=diff.log.bitcoin, main = "Mcleod-Li test for checking changing variance on bitcoin data")

# ARIMA Model Selection

par(mfrow=c(1,2))
acf(diff.log.bitcoin)
pacf(diff.log.bitcoin)

# arch effect visible in the both acf and pacf plots

eacf(diff.log.bitcoin)

# ARIMA(2,1,2), ARIMA(1,1,2), ARIMA(1,1,1)

res = armasubsets(y=diff.log.bitcoin,nar=7,nma=7,y.name='test',ar.method='ols')
plot(res)

# ARIMA(5,1,5),ARIMA(6,1,6) 

# Overall;  ARIMA(2,1,2), ARIMA(1,1,2), ARIMA(1,1,1), ARIMA(5,1,5),ARIMA(6,1,6)

# ARIMA(1,1,1) # not significant 

model_111_css = arima(log(Bitcoin_ts),order=c(1,1,1),method='CSS') 
coeftest(model_111_css) 

model_111_ml = arima(log(Bitcoin_ts),order=c(1,1,1),method='ML') 
coeftest(model_111_ml) 

residual.analysis(model_111_ml, std = TRUE,start = 1)


# ARIMA(1,1,2) # Not significant 

model_112_css = arima(log(Bitcoin_ts),order=c(1,1,2),method='CSS') 
coeftest(model_112_css) 

model_112_ml = arima(log(Bitcoin_ts),order=c(1,1,2),method='ML') 
coeftest(model_112_ml) 

residual.analysis(model_211_ml, std = TRUE,start = 1)

# ARIMA(2,1,2) #All significant

model_212_css = arima(log(Bitcoin_ts),order=c(2,1,2),method='CSS') 
coeftest(model_212_css) 

model_212_ml = arima(log(Bitcoin_ts),order=c(2,1,2),method='ML') 
coeftest(model_212_ml) 

residual.analysis(model_212_ml, std = TRUE,start = 1)

# ARIMA(5,1,5) # not significant

model_515_css = arima(log(Bitcoin_ts),order=c(5,1,5),method='CSS') 
coeftest(model_515_css) 

model_515_ml = arima(log(Bitcoin_ts),order=c(5,1,5),method='ML') 
coeftest(model_515_ml) 

residual.analysis(model_015_ml, std = TRUE,start = 1)

# ARIMA(6,1,6)  #not significant

model_616_css = arima(log(Bitcoin_ts),order=c(6,1,6),method='CSS') 
coeftest(model_616_css) 

model_616_ml = arima(log(Bitcoin_ts),order=c(6,1,6),method='ML') 
coeftest(model_616_ml) 

residual.analysis(model_616_ml, std = TRUE,start = 1)

# Comparing the AIC and BIC results

sort.score(AIC(model_111_ml,model_112_ml,model_212_ml,model_515_ml, model_616_ml), score="aic")
sort.score(BIC(model_111_ml,model_112_ml,model_212_ml,model_515_ml, model_616_ml), score="bic")

#Underfitting check of ARIMA(6,1,6)

# ARIMA(5,1,6) 

model_516_css = arima(log(Bitcoin_ts),order=c(5,1,6),method='CSS') 
coeftest(model_516_css) 

model_516_ml = arima(log(Bitcoin_ts),order=c(5,1,6),method='ML') 
coeftest(model_516_ml) 

residual.analysis(model_516_ml, std = TRUE,start = 1)

# ARIMA(6,1,5) 

model_615_css = arima(log(Bitcoin_ts),order=c(6,1,5),method='CSS') 
coeftest(model_615_css) 

model_615_ml = arima(log(Bitcoin_ts),order=c(6,1,5),method='ML') 
coeftest(model_615_ml) 

residual.analysis(model_615_ml, std = TRUE,start = 1)

# Residual Analysis of ARIMA(6,1,6)

m616_residuals = model_616_ml$residuals

# Absolute Value and Square Root Transformations

abs.res = abs(m616_residuals)
sq.res = m616_residuals^2

# Absolute Value Part
acf(abs.res, ci.type="ma",main="The sample ACF plot for absolute residual series")
pacf(abs.res, main="The sample PACF plot for absolute residual series")
eacf(abs.res)

#Square Root Part
acf(sq.res, ci.type="ma",main="The sample ACF plot for absolute residual series")
pacf(sq.res, main="The sample PACF plot for absolute residual series")
eacf(sq.res)


# From the EACF of absolute residuals, we can identify  ARMA(1,1),ARMA(1,2),ARMA(2,2) models for absolute residual series. 
# These models correspond to parameter settings of [max(1,1),1], [max(1,2),1], and [max(2,2),2]. So the corresponding 
# tentative GARCH models are GARCH(1,1),GARCH(2,1),GARCH(2,2).


# From the EACF of squared residuals, we can identify  ARMA(2,3),ARMA(2,4),ARMA(3,4) models for absolute residual series. 
# These models correspond to parameter settings of [max(2,3),2], [max(2,4),2], and [max(3,4),3]. So the corresponding 
# tentative GARCH models are GARCH(3,2),GARCH(4,2), GARCH(4,3).

# ARMA + GARCH
# ARMA(6,6) and GARCH (1,1)

model1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                   mean.model = list(armaOrder = c(6, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.66_11<-ugarchfit(spec = model1, data = diff.log.bitcoin, out.sample = 100)
m.66_11  
plot(m.66_11)

# ARMA(6,6) and GARCH (2,1)

model2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 1)), 
                   mean.model = list(armaOrder = c(6, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.66_21<-ugarchfit(spec = model2, data = diff.log.bitcoin, out.sample = 100)
m.66_21  
plot(m.66_21)

# ARMA(6,6) and GARCH (2,2)

model3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), 
                   mean.model = list(armaOrder = c(6, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.66_22<-ugarchfit(spec = model3, data = diff.log.bitcoin, out.sample = 100)
m.66_22  

plot(m.66_22)

# ARMA(6,6) and GARCH (4,3)

model4<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 3)), 
                   mean.model = list(armaOrder = c(6, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.66_43<-ugarchfit(spec = model4, data = diff.log.bitcoin, out.sample = 100)
m.66_43  
plot(m.66_43)

# ARMA(6,6) and GARCH (4,2)

model5<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 2)), 
                   mean.model = list(armaOrder = c(6,6), include.mean = FALSE), 
                   distribution.model = "norm")
m.66_42<-ugarchfit(spec = model5, data = diff.log.bitcoin, out.sample = 100)
m.66_42  
plot(m.66_42)

# ARMA(6,6) and GARCH (4,3)

model5<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 3)), 
                   mean.model = list(armaOrder = c(6,6), include.mean = FALSE), 
                   distribution.model = "norm")
m.66_43<-ugarchfit(spec = model5, data = diff.log.bitcoin, out.sample = 100)
m.66_43  
plot(m.66_43)

# Forecasting

forc.66_21 = ugarchforecast(m.66_21, data = diff.log.bitcoin, n.ahead = 10, n.roll = 10)
plot(forc.66_21, which = "all")
forc.66_21

# Reversing data to original 
log.data.diff1.back = diffinv(diff.log.bitcoin, xi = log(Bitcoin_ts)[1])
log.data.diff1.back = exp(log.data.diff1.back)


spec <- ugarchspec(variance.model = list(model ="sGARCH", garchOrder = c(2,1), submodel =NULL,external.regressors =NULL, variance.targeting =FALSE),mean.model= list(armaOrder = c(6,6),external.regressors =NULL,distribution.model ="std",start.pars = list(),fixed.pars = list()))
model.AR_GARCH <- ugarchfit(spec = spec, data = Bitcoin_ts,solver.control = list(trace=0))
model.AR_GARCH
forc <- ugarchforecast(model.AR_GARCH,n.ahead=10,data=Bitcoin_ts)
plot(forc)
forc
a = forc@forecast$seriesFor[,ncol(forc@forecast$seriesFor)]
a


# MASE For Model Fitting

which(is.na(Bitcoin$Close))

Bitcoin$Close[1773] <- 11573.30

Bitcoin$Close[1774] <- 10779.90

Bitcoin$Close[1775] <- 9965.57

Bitcoin$Close[1776] <- 9395.01

Bitcoin$Close[1777] <- 9337.55

Bitcoin$Close[1778] <- 8866.00

Bitcoin$Close[1779] <- 9578.63

Bitcoin$Close[1780] <- 9205.12

Bitcoin$Close[1781] <- 9194.85

Bitcoin$Close[1782] <- 8269.81

which(is.na(Bitcoin$Close))


model_616_ml <- ugarchspec()
m.fit <- ugarchfit(spec = model_616_ml, data = Bitcoin$Close)
fitted.values = fitted(m.fit)

MASE(Bitcoin$Close,fitted.values)

# MASE For Forecast
Bitcoin_Prices_Forecasts <- read_excel("C:/Users/Syed Hassan Afsar/Downloads/RMIT/2nd Semester/Time Series Analysis/Final Project/Bitcoin_Prices_Forecasts.xlsx")
MASE(Bitcoin_Prices_Forecasts$`Closing price`,as.numeric(a))
