library(dynlm)
library(stargazer)
library(pdfetch) 
library(xts)
library(urca)
library(zoo)
library(ggplot2)
library(tidyr)
library(dplyr)
library(strucchange)
library(ARDL)
library(lmtest)
library(car)
library(sandwich)

# FETCHING
debtdata <- pdfetch_FRED("HDTGPDUSQ163N")
unempdata <- pdfetch_FRED("LRUN64TTUSQ156S")

colnames(debtdata) [1] <- "debt"
colnames(unempdata) [1] <- "unemp"

# Make into dataframe + DECOMPOSE
hdf = data.frame(debtdata) 
hdf$t<-c(1:74)
hdf$date <- as.Date(index(debtdata))

# Preliminary visualisation
h_ts <- ts (debtdata, freq = 4)
stl = stl(h_ts [,1], s.window="periodic")
plot(stl) 
hdcomp<-decompose(h_ts)
plot(hdcomp)

# SEASONALITY
hdf_SA <- hdf$debt - hdcomp$seasonal
hdf$hSA <- hdf_SA

# UNEMP data
udf = data.frame(unempdata)
udf$t<-c(1:215)
udf$date <- as.Date(index(unempdata))

# ggplots
plot_debt <- ggplot(hdf) + geom_line(aes(x = date, y = debt,)) + ylim(40,120) +
  labs(title = "Household debt to GDP ratio, USA, quarterly (%)", x = "Year", y = "HH Debt to GDP") + 
  theme_minimal()

plot_unemp <- ggplot(udf) + geom_line(aes(x = date, y = unemp)) + ylim(0,15) +
  labs(title = "Unemployment rate, USA, quarterly (%)", x = "Year", y = "Unemployment rate, %") + 
  theme_minimal()

plot_debt
plot_unemp

# create combined df left joined on household debt
hudf <- left_join(hdf, udf, by = c("date"), unmatched = "drop")

plot_unemp2 <- ggplot(hudf) + geom_line(aes(x = date, y = unemp)) + ylim(0,15) +
  labs(title = "Unemployment rate, USA
       2005 - 2023, quarterly (%)", x = "Year", y = "Unemployment rate, %") + 
  theme_minimal()

plot_unemp2

# TRENDS
# checking for trends in hdf
h_linear_model <- lm(data = hdf, hSA ~ t)
h_quad_model <- lm(data = hdf, hSA ~ t + I(t^2))
h_cubic_model <- lm(data = hdf, hSA ~ t + I(t^2) + I(t^3))

summary(h_linear_model)
summary(h_quad_model)
summary(h_cubic_model)

# finding which model is best
AIC(h_linear_model, h_quad_model, h_cubic_model)
BIC(h_linear_model, h_quad_model, h_cubic_model)

# plotting the trends hdf
plot(hdf$date, hdf$hSA, main="Time Series with Fitted Trends")
lines(hdf$date, predict(h_linear_model), col="blue")
lines(hdf$date, predict(h_quad_model), col="red")
lines(hdf$date, predict(h_cubic_model), col="green")
legend("topright", legend=c("Linear", "Quadratic", "Cubic"), col=c("blue", "red", "green"), lty=1)

# checking for trends in udf
u_linear_model <- lm(data = udf, unemp ~ t)
u_quad_model <- lm(data = udf, unemp ~ t + I(t^2))
u_cubic_model <- lm(data = udf, unemp ~ t + I(t^2) + I(t^3))

summary(u_linear_model)
summary(u_quad_model)
summary(u_cubic_model)

# finding which model is best
AIC(u_linear_model, u_quad_model, u_cubic_model)
BIC(u_linear_model, u_quad_model, u_cubic_model)

#plotting the trends udf
plot(udf$date, udf$unemp, main="Time Series with Fitted Trends")
lines(udf$date, predict(u_linear_model), col="blue")
lines(udf$date, predict(u_quad_model), col="red")
lines(udf$date, predict(u_cubic_model), col="green")
legend("topright", legend=c("Linear", "Quadratic", "Cubic"), col=c("blue", "red", "green"), lty=1)

# obtaining detrended residual data
# debt (hdf)
h_quad_detrended <- residuals(h_quad_model)

# unemp (udf)
u_linear_detrended <- residuals(u_linear_model)

# both residual plots
plot(h_quad_detrended)
plot(u_linear_detrended)

#testing resids for unit root
detrended_h <- ur.df(h_quad_detrended, type = "none")
summary(detrended_h)

detrended_u <- ur.df(u_linear_detrended, type = "none")
summary(detrended_u)

# STATIONARITY
adf1 <- ur.df(hdf$hSA, type="none")
summary(adf1)

adf2 <- ur.df(udf$unemp, type="none")
summary(adf2)

# Addressing stationarity in hSA and unemp
# and testing again for unit roots
hdf$diff <- c(NA,diff(hdf$hSA))
adf3 <- ur.df(na.omit(hdf$diff), type = "none")
summary(adf3)

udf$diff <- c(NA,diff(udf$unemp))
adf4 <- ur.df(na.omit(udf$diff), type = "none")
summary(adf4)

# STRUCTURAL BREAKS in hh debt data
# for trending data
bpdebt <- breakpoints(hdf$hSA ~ 1)
summary(bpdebt)
bpdebt$breakpoints
breakpoints(bpdebt, breaks=4)

# for residuals -> no trends
# linear residuals
bp_h_lr <- breakpoints(h_linear_detrended ~ 1)
breakpoints(bp_h_lr, breaks=3)
summary(bp_h_lr)

# quadratic resids
bp_h_qr <-  breakpoints(h_quad_detrended ~ 1)
breakpoints(bp_h_qr, breaks=3)
summary(bp_h_qr)

# cubic resids
bp_h_cr <- breakpoints(h_cubic_detrended ~ 1)
breakpoints(bp_h_cr, breaks=3)
summary(bp_h_cr)

# Plotted breakpoints
bp <- as.Date(c("2007-09-30", "2011-06-30", "2019-12-31"))

plot_debt_bp <- ggplot(hdf) + geom_line(aes(x = date, y = hSA)) + 
  ylim(40,120) + 
  geom_vline(xintercept = as.numeric(bp), linetype = "dashed", color = "red") +
  labs(title = "Household debt to GDP ratio, USA, quarterly", x = "Year", y = "HH Debt to GDP") + 
  theme_minimal()

plot_debt_bp

# dummies
hudf$s2007 <- ifelse(hudf$t.x >= 11 & hudf$t.x < 26, 1, 0)  
hudf$j2011 <- ifelse(hudf$t.x >= 26 & hudf$t.x < 60, 1, 0)
hudf$d2019 <- ifelse(hudf$t.x >= 60, 1, 0)

# ARDL
ts_hudf <- ts(hudf, freq = 4)

# model given in guidance (equation 1)
ardloriginal <- dynlm(hSA~L(hSA) + unemp, data = ts_hudf)
ardl_textm <- stargazer(ardloriginal, type = "text")
dwtest(ardloriginal)

# adding SB dummies
ardldummies <- dynlm(hSA~L(hSA) + unemp + s2007 + j2011 + d2019, data = ts_hudf)
ardldummies_textm <- stargazer(ardldummies, type = "text")

# testing model with SB dummies for optimal amount of lags with BIC
ardldummies_bic<-auto_ardl(hSA~ unemp + s2007 + j2011 + d2019, selection="BIC", data = ts_hudf, max_order = 5)
ardldummies_bic

# models after testing for optimal lags
ardl2 <- ardl(hSA ~ unemp + s2007 + j2011 + d2019, data = hudf, order=c(6,7,7,7,7))
summary(ardl2)

# Serial Correlation
dwtest(ardl2)
bgtest(ardl2, type = "F")

# Heteroskedasticity
bptest(ardl2)

#FINAL
finalmodel <- dynlm(hSA ~ L(hSA) + L(hSA, 2) + L(hSA, 3) + L(hSA, 4) + L(hSA, 5) + L(hSA, 6) + 
                      unemp + L(unemp) + L(unemp, 2) + L(unemp, 3) +L(unemp, 4) + L(unemp, 5) + L(unemp, 6) + L(unemp, 7) +
                      s2007 + L(s2007) + L(s2007, 2) + L(s2007, 3) + L(s2007, 4) + L(s2007, 5) + L(s2007, 6) + L(s2007, 7) +
                      j2011 + L(j2011) + L(j2011, 2) + L(j2011, 3) + L(j2011, 4) + L(j2011, 5) + L(j2011, 6) + L(j2011, 7) +
                      d2019 + L(d2019) + L(d2019, 2) + L(d2019, 3) + L(d2019, 4) + L(d2019, 5) + L(d2019, 6) + L(d2019, 7)
                    , data = ts_hudf)

textmodel <- stargazer(finalmodel, type = "text")