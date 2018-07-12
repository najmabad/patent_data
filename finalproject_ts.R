## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
#DATA SETUP
patent_data<-read.table('patent_counts_by_semester.csv',header=T, sep='\t')
patent_data$time<-as.Date(patent_data$time)

#Patent H01
ipc_H01<-dplyr::filter(patent_data,ipc=="H01")
ts_H01<-ts(ipc_H01[,"count"],frequency=2,start=c(1852,2))
plot(ts_H01, main="PATSTAT time series - IPC code H01 (1852-2010)", ylab="Volume")
H01_subset<-dplyr::filter(patent_data, ipc=="H01", 
                          time>="1902-07-01", 
                          time<="2010-01-01")
ts_H01_subset<-ts(H01_subset[,"count"],frequency=2,start=c(1902,2))
plot(ts_H01_subset, main="PATSTAT time series - IPC code H01 (1902-2010)", ylab="Volume")

#Patent G06
ipc_G06<-dplyr::filter(patent_data,ipc=="G06")
ts_G06<-ts(ipc_G06[,"count"],frequency=2,start=c(1852,2))
plot(ts_G06,main="PATSTAT time series - IPC code G06 (1852-2010)", ylab="Volume")
G06_subset<-dplyr::filter(patent_data, ipc=="G06", 
                          time>="1902-07-01", 
                          time<="2010-01-01")
ts_G06_subset<-ts(G06_subset[,"count"],frequency=2,start=c(1902,2))
plot(ts_G06_subset, main="PATSTAT time series - IPC code G06 (1902-2010)", ylab="Volume")

## ----echo=FALSE,fig.lp="fig:1", fig.cap = "PATSTAT time series - IPC code H01 and G06 (1902-2010)", fig.pos="H"----

#plots
par(mfrow=c(2,2))
plot(ts_H01_subset, main="Patents H01", ylim=c(0, 113019), yaxt='n', ylab='Volume')
axis(side = 2, labels=c("0", "50000", "100000"), at = c(0, 50000, 100000))

plot(log(ts_H01_subset), main="Log of Patents H01", yaxt='n', ylab='Volume')
axis(side = 2, at = c(5, 6, 7, 8, 9, 10, 11, 12))

plot(ts_G06_subset, main="Patents G06", ylim=c(0, 113019), yaxt='n', ylab='Volume')
axis(side = 2, labels=c("0", "50000", "100000"), at = c(0, 50000, 100000))

plot(log(ts_G06_subset), main="Log of Patents G06", yaxt='n', ylab='Volume')
axis(side = 2, at = c(5, 6, 7, 8, 9, 10, 11, 12))


## ----echo=FALSE, results='hide', fig.lp="fig:2", fig.cap = "Simple Exponetial smoothing obtained through the Holt-Winters algorithm", fig.pos="H"----
#EXPONENTIAL SMOOTHING

HW_H01<-HoltWinters(ts_H01_subset,beta=FALSE,gamma=FALSE)
HW_H01
plot(HW_H01, main="Simple Exponetial Smoothing",ylim=c(0, 113019), yaxt='n', ylab='Volume')
axis(side = 1, at = c(1910, 1930, 1950, 1970, 1990, 2010))
axis(side = 2, labels=c("0", "50000", "100000"), at = c(0, 50000, 100000))

## ----echo=FALSE, results='hide', fig.keep='none', warning=FALSE----------
#MODEL SPECIFICATION 

#acf and pacf of the original series
acf(ts_H01_subset, main="ACF of Patent H01",lag.max=100)
pacf(ts_H01_subset, main="PACF of Patent H01",lag.max=100)
adf.test(ts_H01_subset)

## ----echo=FALSE, results='hide', fig.keep='none', warning=FALSE----------

#log and differenciated series
log_ts_H01_subset<-log(ts_H01_subset)
diff_log_ts_H01_subset<-diff(log_ts_H01_subset)
plot(diff_log_ts_H01_subset)
mean(diff_log_ts_H01_subset)
adf.test(diff_log_ts_H01_subset)

## ----echo=FALSE, fig.lp="fig:3", fig.cap = "Logged and differentiated data, ACF and PACF - H01 subset", fig.pos="H", warning=FALSE----

#acf and pacf of the transformed series
par(mfrow=c(3,1))
plot(diff_log_ts_H01_subset, main="Patent H01, first difference and log transformation", ylab="Volume" )
acf(diff_log_ts_H01_subset, main="ACF of trasformed Patent H01",lag.max=40)
pacf(diff_log_ts_H01_subset, main="PACF of trasformed Patent H01",lag.max=40)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------

#models
mod1<-arima(log_ts_H01_subset,order= c(1,1,0))
mod2<-arima(log_ts_H01_subset,c(2,1,0))
mod3<-auto.arima(log_ts_H01_subset,seasonal = FALSE)
auto.arima(log_ts_H01_subset)
auto.arima(log_ts_H01_subset,seasonal = TRUE)

summary(mod1)
summary(mod2)
summary(mod3)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------

#MODEL CHECKING
#residual diagnostics
tsdiag(mod1)
tsdiag(mod2)
tsdiag(mod3)

#normality of the residuals
par(mfrow = c(3,1)) 
qqnorm(mod1$residuals, main='Model 1 - Normal Q-Q Plot');qqline(mod1$residuals);
qqnorm(mod2$residuals, main='Model 2 - Normal Q-Q Plot');qqline(mod2$residuals);
qqnorm(mod3$residuals, main='Model 3 - Normal Q-Q Plot');qqline(mod3$residuals);

shapiro.test(mod1$residuals)
shapiro.test(mod2$residuals)
shapiro.test(mod3$residuals)

#MSFE and MAE
mse_mod1 = mean((mod1$residuals) ^ 2)
mse_mod2 = mean((mod2$residuals) ^ 2)
mse_mod3 = mean((mod3$residuals) ^ 2)

mae_mod1 = mean(abs(mod1$residuals))
mae_mod2 = mean(abs(mod2$residuals))
mae_mod3 = mean(abs(mod3$residuals))

'Model 1 (MSE/MAE):'; mse_mod1; mae_mod1
'Model 2 (MSE/MAE):'; mse_mod2; mae_mod2
'Model 3 (MSE/MAE):'; mse_mod3; mae_mod3

#information criteria
AIC(mod1,mod2,mod3)
BIC(mod1,mod2,mod3)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
#ALTERNATIVE SPECIFICATION
#breakpoints
bp.logpat <- breakpoints(ts_H01_subset ~ 1)
summary(bp.logpat)
breakdates(bp.logpat)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
#time series(1971-2010)
H01_1971<-dplyr::filter(patent_data, ipc=="H01", 
                        time>="1971-07-01",
                        time<="2010-01-01")

ts_H01_1971<-ts(H01_1971[,"count"],frequency=2,start=c(1971,1))
plot(ts_H01_1971)
adf.test(ts_H01_1971)

diff_ts_H01_1971<-diff((ts_H01_1971))
plot(diff_ts_H01_1971)
adf.test(diff_ts_H01_1971)

log_diff_ts_H01_1971<-diff(log(ts_H01_1971))
plot(log_diff_ts_H01_1971)
adf.test(log_diff_ts_H01_1971)
#model specification
acf(log_diff_ts_H01_1971, lag.max=50)
pacf(log_diff_ts_H01_1971, lag.max=50)
mod4<-arima(log(ts_H01_1971),order=c(2,1,0))

#residuals diagnostic
tsdiag(mod4)
qqnorm(mod4$residuals, main='Model 4 - Normal Q-Q Plot');qqline(mod4$residuals);
shapiro.test(mod4$residuals)

mse_mod4 = mean((mod4$residuals) ^ 2)
mae_mod4 = mean(abs(mod4$residuals))
'Model 4 (MSE/MAE):';mse_mod4;mae_mod4
AIC(mod4)
BIC(mod4)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------

#FORECAST

par(mfrow=c(1,1))
forecast_mod1 <- forecast(mod1,h=2,level=95)
forecast_mod2 <- forecast(mod2,h=2,level=95)
forecast_mod3 <- forecast(mod3,h=2,level=95)
forecast_mod4 <- forecast(mod4,h=2,level=95)
forecast_mod1
forecast_mod2
forecast_mod3
forecast_mod4

## ----echo=FALSE, results='hide', fig.lp="fig:4", fig.cap = "Forecast for the two semesters of 2010 using the H01 subset (1971-2010)", fig.pos="H"----
plot(forecast_mod4,main="Forecast for 2010", ylab = "Volume")

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------

#MLEs of the unknown parameters

log_ts_H01_subset<-log(ts_H01_subset)

build<-function(par){dlmModPoly(order=2, dV=(par[1]),dW=c(0,par[3]))}
fit<-dlmMLE(log_ts_H01_subset, par=c(1,0,1), build,lower=c(10e-7,0,0),hessian=T)
show(fit)
fit$convergence
unlist(build(fit$par)[c("V","W")])
V_H<-unlist(build(fit$par)["V"])
W_H<-unlist(build(fit$par)["W"])
W4_H<-W_H[4]

#DLM
coef<-fit$par
dlm_1<-build(coef)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------

# SMOOTHING THE DATA

smoothed_H01 <- dlmSmooth(log_ts_H01_subset, dlm_1)
variances_smoothed_H01<- dlmSvd2var(smoothed_H01$U.S, smoothed_H01$D.S)

#smoothed estimates of the level state, plot and 95% credible intervals
smoothed_H01_level <- smoothed_H01$s[,1]
sd_smoothed_H01_level  <- sapply(variances_smoothed_H01, function(x) sqrt(x[1,1]))
ci_slH01_upper <- smoothed_H01_level + qnorm(0.975, sd=sd_smoothed_H01_level)
ci_slH01_lower <- smoothed_H01_level - qnorm(0.975, sd=sd_smoothed_H01_level)

plot(log_ts_H01_subset,
     xlab = "", ylab = "Level",
     main="Smoothed estimates of the level with 95% credible intervals")   
lines(dropFirst(smoothed_H01_level), col="red")
lines(dropFirst(ci_slH01_upper), lty=2, col= "green")
lines(dropFirst(ci_slH01_lower), lty=2, col= "green" )
legend(
  "bottomright",
  legend=c("data","smoothed level", "95% credible interval"),
  col=c("black", "red","green"),
  lty=c(1,1,2),
  cex=0.5
)

#smoothed estimates of the growth state, plot and 95% credible intervals
smoothed_H01_growth <- smoothed_H01$s[,2]
sd_smoothed_H01_growth <- sapply(variances_smoothed_H01, function(x) sqrt(x[2,2]))
ci_sgH01_upper <- smoothed_H01_growth + qnorm(0.975, sd=sd_smoothed_H01_growth)
ci_sgH01_lower <- smoothed_H01_growth - qnorm(0.975, sd=sd_smoothed_H01_growth)

plot(smoothed_H01_growth, main="Smoothed estimates of the growth rate with 95% credible intervals")
lines(dropFirst(ci_sgH01_upper), lty=2, col= "green")
lines(dropFirst(ci_sgH01_lower), lty=2, col= "green")
legend(
  "topright",
  legend=c("smoothed growth rate", "95% credible interval"),
  cex=0.5,
  col=c("black", "green"),
  lty=c(1,2)
)

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
# FILTERING

# filtered estimates 
filterted_H01 <-dlmFilter(log_ts_H01_subset, dlm_1)
variances_filtered_H01<-dlmSvd2var(filterted_H01$U.C, filterted_H01$D.C)

#filtered estimates of the level state, plot and 95% credible intervals
filtered_H01_level <- filterted_H01$m[,1]
sd_filtered_H01_level  <- sapply(variances_filtered_H01, function(x) sqrt(x[1,1]))
ci_flH01_upper <- filtered_H01_level + qnorm(0.975, sd=sd_filtered_H01_level)
ci_flH01_lower <- filtered_H01_level - qnorm(0.975, sd=sd_filtered_H01_level)

plot(log_ts_H01_subset,
     xlab = "", ylab = "Level", main="Filtered estimates of the level with 95% credible intervals")   
lines(dropFirst(filtered_H01_level), lty = "longdash",col="red")
lines(dropFirst(ci_flH01_upper), lty=2, col= "green")
lines(dropFirst(ci_flH01_lower), lty=2, col= "green")
legend("bottomright",
       legend=c("data","filtered level", "95% credible interval"),
       cex=0.6,
       col=c("black", "red","green"),
       lty=c(1,1,2))

#filtered estimates of the growth state, plot and 95% credible intervals

filtered_H01_growth <- filterted_H01$m[,2]
sd_filtered_H01_growth <- sapply(variances_filtered_H01, function(x) sqrt(x[2,2]))

ci_fgH01_upper <- filtered_H01_growth + qnorm(0.975, sd=sd_filtered_H01_growth)
ci_fgH01_lower <- filtered_H01_growth - qnorm(0.975, sd=sd_filtered_H01_growth)

plot(filtered_H01_growth,
     xlab = "", ylab = "Level", main="Filtered estimates of the growth with 95% credible intervals")   
lines(ci_fgH01_upper, lty=2, col= "green")
lines(ci_fgH01_lower, lty=2, col= "green")
legend("topright",
       legend=c("filtered growth", "95% credible interval"),
       cex=0.6,
       col=c("black","green"),
       lty=c(1,1,2))

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
# FORECASTING 1-step ahead

forecasted_H01<-filterted_H01$f
variances_forecasted_H01 <- dlmSvd2var(filterted_H01$U.R, filterted_H01$D.R)
sd_forecasted_H01<- sapply(variances_forecasted_H01, function(x) sqrt(x[1,1]))

ci_forecastedH01_upper <-  forecasted_H01 + qnorm(0.975, sd=sd_forecasted_H01)
ci_forecastedH01_lower <- forecasted_H01 - qnorm(0.975, sd=sd_forecasted_H01)

plot(log_ts_H01_subset, main="One-step ahead forecasts with 95% credible intervals")
lines(dropFirst(forecasted_H01), col="red")
lines(dropFirst(ci_forecastedH01_upper),lty=2,col="green")
lines(dropFirst(ci_forecastedH01_lower),lty=2,col="green")
legend("bottomright",
       legend=c("data","one-step ahead forecast", "95% credible intervals"),
       col=c("black", "red","green"),
       lty=c(1,1,2), cex=0.6)


## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
# predictive performance
H01_filt <- dlmFilter(log_ts_H01_subset, dlm_1)
res_H01_filt <- residuals(H01_filt, sd=FALSE)
qqnorm(res_H01_filt)
qqline(res_H01_filt)
shapiro.test(res_H01_filt)
tsdiag(H01_filt)

#MSFE and MAE
forecast_residuals_H01<-(dropFirst(log_ts_H01_subset) - dropFirst(forecasted_H01))
msfe <- mean(forecast_residuals_H01^2)
mafe <- mean(abs(forecast_residuals_H01))
msfe
mafe


## ----echo=FALSE, results='hide', fig.keep='none', warning=FALSE----------
# Model and MLE of the unknown parameters
log_ts_G06_subset<-log(ts_G06_subset)

buildG <- function(par){dlmModPoly(order=2, dV=(par[1]),dW=c(0,par[3]))}
fitG <- dlmMLE(log_ts_G06_subset, par=c(2,0,1), buildG, lower=c(10e-7,0,0), hessian=TRUE)
show(fitG)
fitG$convergence
unlist(buildG(fitG$par)[c("V","W")])
V_G<-unlist(buildG(fitG$par)["V"])
W_G<-unlist(buildG(fitG$par)["W"])
W4_G<-W_G[4]

# smmothed estimates
ModG <- build(fit$par)
smoothtsG <-dlmSmooth(log_ts_G06_subset, ModG)
names(smoothtsG)
s1G <- smoothtsG$s[,1]
R1G <-dlmSvd2var(smoothtsG$U.S, smoothtsG$D.S)
sqrtR1G <- sapply(R1G, function(x) sqrt(x[1,1]))
r1G_upper <- s1G + qnorm(0.975, sd=sqrtR1G)
r1G_lower <- s1G - qnorm(0.975, sd=sqrtR1G)
plot(log_ts_G06_subset,
     xlab = "", ylab = "Level", main="Smoothed estimates of the level with 95% credible interval")   
lines(dropFirst(s1G), lty = "longdash",col="red")
lines(window(r1G_upper,start=start(log_ts_G06_subset)[1]), lty=2, col= "green")
lines(window(r1G_lower,start=start(log_ts_G06_subset)[1]), lty=2, col= "green" )
legend("bottomright",legend=c("data","smoothed level", "95% credible interval"),cex=0.6,
       col=c("black", "red","green"),lty=c(1,1,2))

# growth rate
ModG <- build(fit$par)
smoothtsG <-dlmSmooth(log_ts_G06_subset, ModG)
names(smoothtsG)
s2G <- smoothtsG$s[,2]
RG <-dlmSvd2var(smoothtsG$U.S, smoothtsG$D.S)
sqrtRG <- sapply(RG, function(x) sqrt(x[2,2]))
rG_upper <- s2G + qnorm(0.975, sd=sqrtRG)
rG_lower <- s2G - qnorm(0.975, sd=sqrtRG)
plot(s2G, main="Smoothed estimates of the growth rate with 95% credible interval")
lines(window(rG_upper,start=start(log_ts_G06_subset)[1]), lty=2, col= "green")
lines(window(rG_lower,start=start(log_ts_G06_subset)[1]), lty=2, col= "green")
legend("topright",legend=c("smoothed growth rate", "95% credible interval"),cex=0.6,
       col=c("black", "green"),lty=c(1,2))

# point 2

# filtering the data

# level
ModG <- build(fit$par)
filtertsG <-dlmFilter(log_ts_G06_subset, ModG)
names(filtertsG)
m1G <- filtertsG$m[,1]
CG <-dlmSvd2var(filtertsG$U.C, filtertsG$D.C)
sqrtC1G <- sapply(CG, function(x) sqrt(x[1,1]))
c1G_upper <- m1G + qnorm(0.95, sd=sqrtC1G)
c1G_lower <- m1G - qnorm(0.95, sd=sqrtC1G)
plot(log_ts_G06_subset,
     xlab = "", ylab = "Level", main="Filtered estimates of the level with 95% credible interval")   
lines(dropFirst(m1G), lty = "longdash",col="red")
lines(window(c1G_upper,start=start(log_ts_G06_subset)[1]), lty=2, col= "green")
lines(window(c1G_lower,start=start(log_ts_G06_subset)[1]), lty=2, col= "green" )
legend("bottomright",legend=c("data","filtered level", "95% credible interval"),cex=0.6, col=c("black", "red","green"),lty=c(1,1,2))

# growth rate
ModG <- build(fit$par)
filtertsG <-dlmFilter(log_ts_G06_subset, ModG)
m2G <- filtertsG$m[,2]
C2G <-dlmSvd2var(filtertsG$U.C, filtertsG$D.C)
sqrtC2G <- sapply(C2G, function(x) sqrt(x[2,2]))
c2G_upper <- m2G + qnorm(0.975, sd=sqrtC2G)
c2G_lower <- m2G - qnorm(0.975, sd=sqrtC2G)
plot(m2G)
lines(window(c2G_upper,start=start(ts_G06_subset)[1]), lty=2, col= "green")
lines(window(c2G_lower,start=start(ts_G06_subset)[1]),lty=2, col="green" )

# forecasting 1-step ahead

f1G <-filtertsG$f
listFG <- dlmSvd2var(filtertsG$U.R, filtertsG$D.R)
sqrtFG <- sapply(listFG, function(x) sqrt(x[1,1]))
c2G_upper <- f1G + qnorm(0.975, sd=sqrtFG)
c2G_lower <- f1G - qnorm(0.975, sd=sqrtFG)

plot(dropFirst(log_ts_G06_subset), main="One-step ahead forecasts with 95% credible interval")
lines(dropFirst(f1G), col="red")
lines(dropFirst(c2G_upper),lty=2, col="green" )
lines(dropFirst(c2G_lower),lty=2, col="green" )
legend("bottomright",legend=c("data","one-step-ahead forecast", "95% credible intervals"),cex=0.6, col=c("black", "red","green"),lty=c(1,1,2))

#predictive performance
G06_filt <- dlmFilter(log_ts_G06_subset, ModG)
res_G06_filt <- residuals(G06_filt, sd=FALSE)
qqnorm(res_G06_filt);qqline(res_G06_filt)
shapiro.test(res_G06_filt)
tsdiag(G06_filt)

#MSFE and MAE
forecast_residuals_G06<-(dropFirst(log_ts_G06_subset) - dropFirst(f1G))
msfe <- mean(forecast_residuals_G06^2)
mafe <- mean(abs(forecast_residuals_G06))
msfe
mafe

## ----echo=FALSE, results='hide', fig.keep='none'-------------------------
# Plot H01 & G06
plot(log_ts_H01_subset, main="H01 and G06")
lines(log_ts_G06_subset,col="blue")
legend("bottomright",legend=c("H01","G06"),cex=0.6,
       col=c("black", "blue"),lty=c(1,1))

# Merge the two time series
ts_merged<-cbind(log_ts_H01_subset,log_ts_G06_subset)
unlist(ts_merged)

# SUTSE
sutse_dlm <- function ( parameters ){
  # V matrix
  Vsd <- exp( parameters[1:2])
  Vcorr <- tanh(parameters[6])
  V <- Vsd %o% Vsd
  V[1 ,2] <- V[2 ,1] <- V[1 ,2] * Vcorr
  # W matrix
  Wsd <- exp( parameters[3:4])
  Wcorr <- tanh(parameters[5])
  W <- Wsd %o% Wsd
  W[1 ,2] <- W[2 ,1] <- W[1 ,2] * Wcorr
  mod_basis <- dlmModPoly(2)
  mod_basis$FF <- mod_basis$FF %x% diag(2)
  mod_basis$GG <- mod_basis$GG %x% diag(2)
  mod_basis$m0 <- rep(0 ,4)
  mod_basis$C0 <- 1e7 * diag(4)
  mod_basis$V <- V
  mod_basis$W <- bdiag ( matrix (0, ncol =2, nrow =2) , W)
  return (mod_basis )
}
initial_vals <- c(V_H,V_G,W4_H,W4_G,0.1,0.1)
sutse_mle <- dlmMLE (ts_merged, initial_vals , sutse_dlm); sutse_mle
names(sutse_mle)



# point(2)

# filtering the data
Mod_sutse <- sutse_dlm(sutse_mle$par)
filter_sutse <-dlmFilter(ts_merged, Mod_sutse)
names(filter_sutse)
m1 <- filter_sutse$m[,1]
C1 <-dlmSvd2var(filter_sutse$U.C, filter_sutse$D.C)
sqrtC1 <- sapply(C1, function(x) sqrt(x[1,1]))
c1_upper <- m1 + qnorm(0.95, sd=sqrtC1)
c1_lower <- m1 - qnorm(0.95, sd=sqrtC1)
m2 <- filter_sutse$m[,2]
C2 <- dlmSvd2var(filter_sutse$U.C, filter_sutse$D.C)
sqrtC2 <- sapply(C2, function(x) sqrt(x[2,2]))
c2_upper <- m2 + qnorm(0.95, sd=sqrtC2)
c2_lower <- m2 - qnorm(0.95, sd=sqrtC2)

plot(log_ts_H01_subset,xlab = "Time", ylab = "", main="Filtered estimates (SUTSE) with 95% credible interval")
lines(log_ts_G06_subset,col="blue")
lines(dropFirst(m1), lty = "longdash",col="red")
lines(window(c1_upper,start=start(ts_merged)[1]), lty=2, col= "green")
lines(window(c1_lower,start=start(ts_merged)[1]), lty=2, col= "green" )
lines(dropFirst(m2), lty = "longdash",col="darkred")
lines(window(c2_upper,start=start(ts_merged)[1]), lty=2, col= "darkgreen")
lines(window(c2_lower,start=start(ts_merged)[1]), lty=2, col= "darkgreen" )
legend("bottomright",legend=c("H01","G06","filtered level for H01","filtered level for G06", "95% credible interval for H01","95% credible interval for G06"),cex=0.6,
       col=c("black","blue","red","darkred","green","darkgreen"),lty=c(1,1,1,1,2,2))

# forecasting 1-step ahead
names(filter_sutse)
unlist(filter_sutse$f)
f1 <-filter_sutse$f[,1]
F1 <- dlmSvd2var(filter_sutse$U.R, filter_sutse$D.R)
sqrtF1 <- sapply(F1, function(x) sqrt(x[1,1]))
f2 <-filter_sutse$f[,2]
F2 <- dlmSvd2var(filter_sutse$U.R, filter_sutse$D.R)
sqrtF2 <- sapply(F2, function(x) sqrt(x[2,2]))

plot(log_ts_H01_subset,ylab="", main="One-step ahead forecasts (SUTSE) with 95% credible interval")
lines(log_ts_G06_subset,col="blue")
lines(f1, col="red")

lines(f1 + 2*(sqrtF1),lty=2,col="green")
lines(f1 - 2*(sqrtF1),lty=2,col="green")
lines(f2, col="darkred")
lines(f2 + 2*(sqrtF1),lty=2,col="darkgreen")
lines(f2 - 2*(sqrtF1),lty=2,col="darkgreen")
legend("bottomright",legend=c("H01","G06","one-step ahead forecasts for H01","one-step ahead forecasts for G06", "95% credible interval for H01","95% credible interval for G06"),cex=0.6,
       col=c("black","blue","red","darkred","green","darkgreen"),lty=c(1,1,1,1,2,2))

# smoothing the data (not required)
Mod_sutse <- sutse_dlm(sutse_mle$par)
smoother_sutse <-dlmSmooth(ts_merged, Mod_sutse)
names(smoother_sutse)
s1 <- smoother_sutse$s[,1]
R1 <-dlmSvd2var(smoother_sutse$U.S, smoother_sutse$D.S)
sqrtR1 <- sapply(R1, function(x) sqrt(x[1,1]))
r1_upper <- s1 + qnorm(0.95, sd=sqrtR1)
r1_lower <- s1 - qnorm(0.95, sd=sqrtR1)
s2 <- smoother_sutse$s[,2]
R2 <-dlmSvd2var(smoother_sutse$U.S, smoother_sutse$D.S)
sqrtR2 <- sapply(R2, function(x) sqrt(x[2,2]))
r2_upper <- s2 + qnorm(0.95, sd=sqrtR2)
r2_lower <- s2 - qnorm(0.95, sd=sqrtR2)

plot(log_ts_H01_subset,xlab = "Time", ylab = "", main="Smoothed estimates (SUTSE) with 95% credible interval")
lines(log_ts_G06_subset,col="blue")
lines(dropFirst(s1), lty = "longdash",col="red")
lines(window(r1_upper,start=start(ts_merged)[1]), lty=2, col= "green")
lines(window(r1_lower,start=start(ts_merged)[1]), lty=2, col= "green" )
lines(dropFirst(s2), lty = "longdash",col="darkred")
lines(window(r2_upper,start=start(ts_merged)[1]), lty=2, col= "darkgreen")
lines(window(r2_lower,start=start(ts_merged)[1]), lty=2, col= "darkgreen" )
legend("bottomright",legend=c("H01","G06","smoothed level for H01","smoothed level for G06", "95% credible interval for H01","95% credible interval for G06"),cex=0.6,
       col=c("black","blue","red","darkred","green","darkgreen"),lty=c(1,1,1,1,2,2))

