## Removing previous object in the R-session
rm(list = ls(all= TRUE))

## loding extra libraries I need for the following analyses
library(glmmADMB)
library(splines)
library(splines)  
library(lme4) 
library(mgcv)

source("~/uio/phenology/phenology.scripts/data.R", echo=T)

phen$transectnr = as.factor(phen$transectnr)
phen$trans = as.factor(phen$trans)
phen$year.x = as.factor(phen$year.x)

#### Exploring the data ####
plot(phen[phen$year.x==2011, ]$Hum~phen[phen$year.x==2011, ]$date.x)
plot(phen[phen$year.x==2011, ]$Hum~phen[phen$year.x==2011, ]$date.x)
plot(phen[phen$year.x==2011, ]$w,phen[phen$year.x==2011, ]$date.x)
plot(phen[phen$year.x==2011, ]$temp ~ns(phen[phen$year.x==2011, ]$date.x, df=2))


plot(phen[phen$year.x==2012, ]$Hum~phen[phen$year.x==2012, ]$date.x)
plot(phen[phen$year.x==2012, ]$Hum~phen[phen$year.x==2012, ]$date.x)
plot(phen[phen$year.x==2012, ]$w,phen[phen$year.x==2012, ]$date.x)


plot(phen[phen$Temp<6,]$tot, phen[phen$Temp<6,]$date)

plot(phen$Hum ~ phen$w)
hist(phen$Hum)
hist(phen$VPD)
nrow(phen[phen$VPD<6,])
nrow(phen[phen$Hum>60,])
nrow(phen[phen$w=="R",])


nrow(phen)
nrow(phen2)
(nrow(phen2[phen2$w=="R",]))/(nrow(phen2))
(nrow(phen[phen$w=="R",]))/(nrow(phen))




phen$datonr2 = phen$date.x
phen$datonr = phen$datonr2
phen[phen$date.x> 370,]$datonr = phen[phen$date.x> 370,]$datonr2-365

#### Checking for trends in proportions between instar stages and covariates, and to see if they are correlations between stages #### 
cor.test((phen$tick_adult_M + phen$tick_adult_F), phen$tick_nymph, method = "spearman") # rho = 0.2748885, p-value < 2.2e-16
cor.test((phen$tick_adult_M + phen$tick_adult_F), phen$tick_nymph, method = "kendall")  # tau = 0.2638104, p-value < 2.2e-16
## We see that both correlations are very weak, but they are correlated - indicating that many nymphs give many adults
## Lets try to see if number of adults can be an effect of number of ticks i a glm-system using glmmADMB. I try this both with and without random effects and choose the models best AIC, and compare results

AIC(adult.mod <- glmmadmb(I(tick_adult_F + tick_adult_M) ~ tot + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) # AIC = 1434.684
AIC(adult.mod.no_rand <- glmmadmb(I(tick_adult_F + tick_adult_M) ~ tot, family="nbinom", data= phen, zeroInflation =T)) # AIC = 1471.77
summary(adult.mod)  ## The best model is the one with random effects
=T)))
# Summary of the best adult model yields a strongly significant effect of total counts of ticks.

## Continue to see if there are any trends in proportions between adult ticks and the covariates. I will report the models with the best AIC- values. These models are optimized using REML - which is the method recommended by Zuur et al. 2009 for comparison of models with different random effects. 

AIC(modbin1 <- gam(cbind(tick_nymph,tot) ~ Temp + s(uniktrans, bs ="re"),  data=phen, family ="binomial", method="REML"))
AIC(modbin1 <- gam(cbind(tick_nymph,tot) ~ Temp,  data=phen, family ="binomial", method="REML"))

AIC(modbin2 <- gam(cbind(tick_nymph,tot) ~ ns(Temp, 3) + s(uniktrans, bs ="re"),  data=phen, family ="binomial", method="REML"))
AIC(modbin2 <- gam(cbind(tick_nymph, tot) ~ ns(Temp, df = 3),  data=phen, family ="binomial", method="REML")) 

AIC(modbin3 <- gam(cbind(tick_nymph, tot) ~ ns(date.x, df=4) + s(uniktrans, bs ="re"),  data=phen, family ="binomial", method="REML"))
AIC(modbin3 <- gam(cbind(tick_nymph, tot) ~ ns(date.x, df=4),  data=phen, family ="binomial", method="REML"))

AIC(modbin4 <- gam(cbind(tick_nymph, tot) ~ alt.cat + s(uniktrans, bs ="re"),  data=phen, family ="binomial", method="REML"))      
AIC(modbin4 <- gam(cbind(tick_nymph, tot) ~ alt.cat,  data=phen, family ="binomial", method="REML"))         
        
AIC(modbin5 <- gam(cbind(tick_nymph, tot) ~ Hum + s(uniktrans, bs ="re"),  data=phen, family ="binomial", method="REML"))        
AIC(modbin5 <- gam(cbind(tick_nymph, tot) ~ Hum,  data=phen, family ="binomial", method="REML"))     

## p-values ranging between: 9.56 and 0.402, and the models without random effects were alwyas the best in AIC
## Se no point in dividing 


## Starting by ploting the weather effects against tick counts to see if there are any unlinearities
plot(phen$Temp, phen$tot)  ## seems to be a second order polynomial
plot(phen$VPD, phen$tot)   ## seems to be a linear effect
plot(phen$Hum, phen$tot)   ## seems to be a linear effect

## Try correlation tests to see if any of the weather variables are correlated, and find predictors that can not be combined as covariates in the models
cor.test(phen$Temp, phen$VPD, method = "spearman") # rho = 0.7362433, p-value < 2.2e-16 - medium correlation, be careful
cor.test(phen$Temp, phen$Hum, method = "spearman") # rho =-0.5622607, p-value < 2.2e-16 - weak correlation, this is ok 
cor.test(phen$Hum,  phen$VPD, method = "spearman") # rho =-0.9686004, p-value < 2.2e-16 - strong correlation, do not combine
cor.test(phen$Temp, phen$VPD, method = "kendall")  # tau = 0.5490592, p-value < 2.2e-16 - medium correlation, be careful
cor.test(phen$Temp, phen$Hum, method = "kendall")  # tau =-0.3942358, p-value < 2.2e-16 - weak correlation, this is ok
cor.test(phen$Hum,  phen$VPD, method = "kendall")  # tau =-0.850031,  p-value < 2.2e-16 - strong correlation, do not combine

#### Model selection ####
##  I first start the model selection by performing three varieties of the full model, namely the using either humidity or vapor pressure deficit as humidity measure, with both 2.order polynomial, and as linear predictors. Temperature seems to be second order from the plot above, and is therefore only tried as a 2.order when I try to fit the optimal random effect.
## See no biological meaning in 

AIC(M1 <- glmmadmb(tot ~ year.x + poly(Hum,2) * poly(Temp, 2) + (1|Sted/trans), family="nbinom", data= phen, zeroInflation =T))  ## AIC= 5967.38
AIC(M2 <- glmmadmb(tot ~ year.x + poly(Hum,2) * poly(Temp, 2) + (1|Sted/uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5967.38
AIC(M3 <- glmmadmb(tot ~ year.x + poly(Hum,2) * poly(Temp, 2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5965.38 


AIC(M4 <- glmmadmb(tot ~ year.x + poly(VPD,2) * poly(Temp, 2) + (1|Sted/trans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5972.06
AIC(M5 <- glmmadmb(tot ~ year.x + poly(VPD,2) * poly(Temp, 2) + (1|Sted/uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5972.06
AIC(M6 <- glmmadmb(tot ~ year.x + poly(VPD,2) * poly(Temp, 2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5970.06 
## Testing some models for all interactions - to se if fit gets more than 2 AIC-points worse when they are removed
AIC(M7 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + Temp:Hum + I(Temp^2):Hum + I(Hum^2):Temp + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5963.4
AIC(M8 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + Temp:Hum + I(Temp^2):Hum + I(Hum^2):I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5963.5
AIC(M9 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + Temp:Hum + I(Hum^2):Temp + I(Hum^2):I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5963.42
AIC(M10 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + I(Temp^2):Hum + I(Hum^2):Temp + I(Hum^2):I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5963.74
AIC(M11 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + Temp:Hum + I(Temp^2):Hum + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5965.22 
AIC(M12 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + Temp:Hum + I(Hum^2):Temp + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5969.04
AIC(M13 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + I(Temp^2):Hum + I(Hum^2):Temp + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5967.04 
AIC(M14 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + Temp:Hum + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5967.08
AIC(M15 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + I(Temp^2):Hum + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5966.16
AIC(M16 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5965.48

## The model without interactions are within Burnham and Andersons 2 units in AIC - and the model are equally fit - I chose the model with the fewes parameters, that is the one without interactions

## Check for effects of 2.order polynomials
AIC(M17 <- glmmadmb(tot ~ year.x + Hum + I(Hum^2) + Temp + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5982.52
AIC(M18 <- glmmadmb(tot ~ year.x + Hum + Temp + I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5964.9 This is the best model
AIC(M19 <- glmmadmb(tot ~ year.x + Hum + Temp + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5980.86

## M18, The model that includes a second order polynomial function of Temperature is the best. Try to remove other variables
AIC(M20 <- glmmadmb(tot ~ year.x + Hum + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5983.48
AIC(M21 <- glmmadmb(tot ~ year.x + Temp + I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC=  5993.62
AIC(M22 <- glmmadmb(tot ~ Hum + Temp + I(Temp^2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5988.9
## Model M18 is still the best.
## Continue by checking if there are any significant time trends in the residuals:
# First, we must make a time variable that has the same value for each day - by subtracting the days for year == 2012 by 365


datodata = phen$date.x
datodata[phen$year.x == 2012] = datodata[phen$year.x == 2012]-365
min(datodata)
max(datodata)
## All is well here - continue by checking residuals 
residualer = resid(M18)   ## Extracting scaled pearson residuals from the best model
## See if there are any trends in the residuals
summary(lm(residualer ~ datodata * alt.cat, data = phen))  ## linear time
summary(lm(residualer ~ ns(datodata, df = 2) * alt.cat, data = phen))  ## Patterns with a spline function with df=2 is significant
summary(lm(residualer ~ ns(datodata, df = 3) * alt.cat, data = phen))  ## No improvement beyond df=2

## We find strong residuals in the trend. Will now try to see if residual trend is removed when we include the altitude category:
AIC(M23 <- glmmadmb(tot ~ year.x + Hum + poly(Temp,2) + alt.cat * ns(datodata, df=2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5840.
AIC(M24 <- glmmadmb(tot ~ year.x + Hum + poly(Temp,2) + alt.cat * datodata + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5882.66
AIC(M25 <- glmmadmb(tot ~ year.x + Hum + poly(Temp,2) + alt.cat + ns(datodata, df=2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5874.28
AIC(M26 <- glmmadmb(tot ~ year.x + Hum + poly(Temp,2) + alt.cat + datodata + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5880.94

## The model including alt.cat and the interaction with ns(date, df=2) is the best model based on AIC
## Checking for residual trend

residualer = resid(M23)   ## Extracting scaled pearson residuals from the best model
## See if there are any trends in the residuals
summary(lm(residualer ~ datodata * alt.cat, data = phen))  ## linear time
summary(lm(residualer ~ ns(datodata, df = 2) * alt.cat, data = phen))  ## Patterns is lost - all is well

#########################################################
######################################################### 


# Try to remove some more weather variables to see if the time trend removes significance in any variables
AIC(M26 <- glmmadmb(tot ~ year.x + Hum + Temp + alt.cat * ns(datodata, df=2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 5848.12
AIC(M27 <- glmmadmb(tot ~ year.x + poly(Temp,2) + alt.cat * ns(datodata, df=2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 
AIC(M28 <- glmmadmb(tot ~ Hum + poly(Temp,2) + alt.cat * ns(datodata, df=2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 
AIC(M29 <- glmmadmb(tot ~ year.x + Hum + alt.cat * ns(datodata, df=2) + (1|uniktrans), family="nbinom", data= phen, zeroInflation =T)) ## AIC= 
## Humidity must not out, and the new best model is still M23



residualer = resid(M27)   ## Extracting scaled pearson residuals from the best model
## See if there are any trends in the residuals
summary(lm(residualer ~ datodata * alt.cat, data = phen))  ## linear time
summary(lm(residualer ~ ns(datodata, df = 2) * alt.cat, data = phen))  ## Patterns is lost - all is well


#Best model is:
     summary(M27)
## Fill into paper!!

#save(list=ls(all=T), file="29.11climateworkspace.RData")


#### Next step is to predict the model, and plot the trend data that from weather model ####
     
model.best = 
newdat.Mspline <- data.frame(expand.grid(datodata = c(min(datodata):max(datodata)), alt.cat = levels(phen$alt.cat), year.x=c(2011, 2012))) 
mm <- model.matrix(delete.response(terms(model.best)),newdat.Mspline)
newdat.Mspline$tot <-(mm %*% fixef(model.best))

plot(newdat.Mspline[newdat.Mspline$alt.cat == "low",]$tot~newdat.Mspline[newdat.Mspline$alt.cat == "low",]$datodata, type ="n")
#points(newdat.Mspline[newdat.Mspline$alt.cat == "low",]$tot~newdat.Mspline[newdat.Mspline$alt.cat == "low",]$datodata)
lines(newdat.Mspline[newdat.Mspline$alt.cat == "low",]$tot~newdat.Mspline[newdat.Mspline$alt.cat == "low",]$datodata)




model.best = M13A
newdat.M23 <- data.frame(expand.grid(Temp=min(phen$Temp):max(phen$Temp), year.x=c("2011", "2012"), Hum = 0))

intercept(sum(ranef(M13A)$'Sted:trans'[1:3])+sum(ranef(M13A)$'Sted:trans'[7:9]))/6


mm <- model.matrix(delete.response(terms(model.best)),newdat.M23)
newdat.M23$tot <-(mm %*% fixef(model.best))
plot(phen$Temp,phen$tot, main = "Tick counts as a function of temperature", ylab = "Tick counts", xlab = "Tempearture", las=1)
lines(exp(newdat.M23[newdat.M23$year.x=="2011",]$tot+2)~newdat.M23[newdat.M23$year.x=="2011",]$Temp, col = "red")
lines(exp(newdat.M23[newdat.M23$year.x=="2012",]$tot+intercept)~newdat.M23[newdat.M23$year.x=="2012",]$Temp, col = "blue")
                 legend(25, 60, c("2011", "2012"),
                 col = c("red", "blue"),
                        lty=1
)
#curve(exp( 20 + 8.3772*x - 13.6844*(x^2)), add=T, col = "red")



save(list=ls(all=T), file="14.11climateworkspace.RData")


M0.11 = glmmadmb(tot ~ 1 + (1|Sted/trans),family="nbinom", data= phen[phen$year.x==2011,])
M0.12 = glmmadmb(tot ~ 1 + (1|Sted/trans),family="nbinom", data= phen[phen$year.x==2012,])

clim11 = glmmadmb(tot ~  Hum * poly(Temp, 2, raw = T) + (1|Sted/trans), family="nbinom", data= phen[phen$year.x==2011,])
clim12 = glmmadmb(tot ~  Hum * poly(Temp, 2, raw = T) + (1|Sted/trans), family="nbinom", data= phen[phen$year.x==2012,])

library(splines)

spline11 <- glmmadmb(tot ~ ns(date.x, df=5) * alt.cat + (1|Sted/trans), data=phen[phen$year.x==2011,], family="nbinom")
spline12 <- glmmadmb(tot ~ ns(date.x, df=5) * alt.cat + (1|Sted/trans), data=phen[phen$year.x==2012,], family="nbinom")
AIC(spline11)
AIC(spline12)

        
        dim(phen[phen$year.y==2011,]$Temp)
       
## artial regressions using only fixed effects from the models
    ## 2011 
    ## Starting with extracting predicted vectors
    model.best = spline11
    newdat.spline11 <- data.frame(date.x=phen[phen$year.x == 2011, ]$date.x, alt.cat = phen[phen$year.x == 2011, ]$alt.cat)
    mm <- model.matrix(delete.response(terms(model.best)),newdat.spline11)
    newdat.spline11$tot <-(mm %*% fixef(model.best))
    
        model.best = clim11
    newdat.clim11 <- data.frame(Temp=phen[phen$year.x == 2011, ]$Temp, Hum = phen[phen$year.x == 2011, ]$Hum)
    mm <- model.matrix(delete.response(terms(model.best)),newdat.clim11)
    newdat.clim11$tot <-(mm %*% fixef(model.best))
plot(phen$tot~phen$Temp)
points(exp(newdat.clim11$tot+3)~newdat.clim11$Temp, col="red")   


    model.best = spline12
    newdat.spline12 <- data.frame(date.x=phen[phen$year.x == 2012, ]$date.x, alt.cat = phen[phen$year.x == 2012, ]$alt.cat)
    mm <- model.matrix(delete.response(terms(model.best)),newdat.spline12)
    newdat.spline12$tot <-(mm %*% fixef(model.best))
    
    model.best = clim12
w=c(as.factor(phen[phen$year.x == 2012, ]$w))
    w[w==1] <-"O"
    w[w==2] <-"R"
w
newdat.clim12 <- data.frame(Temp=phen[phen$year.x == 2012, ]$Temp, Hum = phen[phen$year.x == 2012, ]$Hum)
    mm <- model.matrix(delete.response(terms(model.best)),newdat.clim12)
    newdat.clim12$tot <-(mm %*% fixef(model.best))
        
summary(comb11<-lm(newdat.spline11$tot ~ newdat.clim11$tot-1))
summary(comb12<-lm(newdat.spline12$tot ~ newdat.clim12$tot-1))
    par(mfrow=c(2,2))  
    plot(comb11)
    plot(comb12)
    



scale.res = function(y,mod)
{
  (y-fitted(mod))/sqrt(fitted(mod)+ (fitted(mod))^2)/mod$alpha
}

(resid(M16))
sum((phen$tot-fitted(M16))/sqrt(fitted(M16)))
  

## Plotting temperatures against the predicted vectors for tick abundance 
## 2011
par(mfrow = c(2,1))

plot(phen[phen$year.x==2011,]$temp~phen[phen$year.x==2011,]$date.x, main="Temperature and abundance trends 2011", xlab="Time", ylab="temperature", las=1, type="n")
points(phen[phen$year.x==2011&phen$alt.cat=="high",]$temp~phen[phen$year.x==2011&phen$alt.cat=="high",]$date.x, col="red", pch = 15)
points(phen[phen$year.x==2011&phen$alt.cat=="low",]$temp~phen[phen$year.x==2011&phen$alt.cat=="low",]$date.x, col = "blue", pch = 25)
curve(6.5+ 0*x, lty=2, add=T)
curve(15.5+ 0*x, lty=2, add=T)
#points(exp(comp.vect11)~phen[phen$year.x==2011,]$date.x)
model.best = spline11
newdat <- expand.grid(date.x=c(117:326), alt.cat=c("high","low"))
mm <- model.matrix(delete.response(terms(model.best)),newdat)
newdat$tot <-(mm %*% fixef(model.best))
    
newdat$tot=exp(newdat$tot)
par(new=TRUE)
plot(newdat$date[newdat$alt.cat == "low"&newdat$date>117&newdat$date<326], newdat$tot[newdat$alt.cat == "low"&newdat$date>117&newdat$date<326],, type = "l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", ylim =c(0,4))
axis(4)
lines(newdat$date[newdat$alt.cat == "high"&newdat$date>156&newdat$date<299], newdat$tot[newdat$alt.cat == "high"&newdat$date>156&newdat$date<299], lwd=2, col ="red")

## 2012
plot(phen[phen$year.x==2012,]$temp~phen[phen$year.x==2012,]$date.x, main="Temperature and abundance trends 2012", xlab="Time", ylab="temperature", las=1, type="n")
points(phen[phen$year.x==2012&phen$alt.cat=="high",]$temp~phen[phen$year.x==2012&phen$alt.cat=="high",]$date.x, col="red", pch = 15)
points(phen[phen$year.x==2012&phen$alt.cat=="low",]$temp~phen[phen$year.x==2012&phen$alt.cat=="low",]$date.x, col = "blue", pch=25)
curve(6.5+ 0*x, lty =2, add=T)
curve(15.5+ 0*x, lty=2, add=T)
#points(exp(comp.vect11)~phen[phen$year.x==2011,]$date.x)
model.best = spline12
newdat <- expand.grid(date.x=c(482:662), alt.cat=c("high","low"))
mm <- model.matrix(delete.response(terms(model.best)),newdat)
newdat$tot <-(mm %*% fixef(model.best))
newdat$tot=exp(newdat$tot)
par(new=TRUE)
plot(newdat$date[newdat$alt.cat == "low"&newdat$date>482&newdat$date<662], newdat$tot[newdat$alt.cat == "low"&newdat$date>482&newdat$date<662],, type = "l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", ylim =c(0,4), xlim=c(482,662))
axis(4)
lines(newdat$date[newdat$alt.cat == "high"&newdat$date>509&newdat$date<649], newdat$tot[newdat$alt.cat == "high"&newdat$date>509&newdat$date<649], lwd=2, col ="red")

legend(620, 4, c("Low - temp", "High - temp", "Low - abund", "High - abund", "Thresholds"),
       text.col = "black",
       cex= 0.75,
       lty = c(NA, NA, 1,1,2), 
       pch = c(25, 15, NA,NA,NA),
       col = c("blue", "red", "blue", "red", "black"),
       bg = "white"                        
)

dev.copy2pdf(device=x11,
             #family=c("/home/lars/arial.afm",
             #"/home/lars/arialbd.afm",
             #"/home/lars/ariali.afm",
             #"/home/lars/arialbi.afm"),
             file = "temp.abund.pdf",
             paper="special",
             width=10,
             height=8,
            # horizontal=FALSE
)

#       
    
    
## Plotting temperature spline model.    
## First, making the plot
## 2012
plot(phen$tot~phen$Temp, main="Temperature and abundance trends 2012", ylab="Numbers of ticks", xlab="temperature", las=1, type="n")
points(phen$tot~phen$Temp, col="red", pch = 15)
#points(exp(comp.vect11)~phen[phen$year.x==2011,]$date.x)
model.best = mspline
newdat <- expand.grid(Temp=(c(0:60)/2))
mm <- model.matrix(delete.response(terms(model.best)),newdat)
newdat$tot <-(mm %*% fixef(model.best))
newdat$tot=exp(newdat$tot+3)
#par(new=TRUE)
#plot(newdat$date[newdat$alt.cat == "low"&newdat$date>482&newdat$date<662], newdat$tot[newdat$alt.cat == "low"&newdat$date>482&newdat$date<662],, type = "l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", ylim =c(0,4), xlim=c(482,662))
#axis(4)
lines(newdat$Temp, newdat$tot, lwd=2, col ="red")
    
    legend(620, 4, c("Low - temp", "High - temp", "Low - abund", "High - abund", "Thresholds"),
           text.col = "black",
           cex= 0.75,
           lty = c(NA, NA, 1,1,2), 
           pch = c(25, 15, NA,NA,NA),
           col = c("blue", "red", "blue", "red", "black"),
           bg = "white"                        
    )    
                                           
### Ymse småting som må sjekkes
       phen[phen$tot>0&phen$Temp<5,30:50]
       sum(phen[phen$Temp<5,]$tot)                                 
       sum(phen[phen$Temp<5.5,]$tot)                                        
       sum(phen[phen$Temp<6,]$tot)                                                                            
       sum(phen[phen$Temp<6.5,]$tick_adult_F)                                                                            
                                           
       sum(phen[phen$Temp<7,]$tot)                  
attach(phen[phen$year==2011,])                                           
       plot((phen$tick_adult_M + phen$tick_adult_F)/phen$tick_nymph~phen$date.x)  
summary(testmod <- glmmadmb(I(phen$tick_adult_M + phen$tick_adult_F) ~phen$tick_nymph,family="nbinom"))
                                           
detach(phen[phen$year==2011,])           
                                           
cor.test((phen$tick_adult_M + phen$tick_adult_F), phen$tick_nymph, method = "kendall")

                                           names(phenuc)
        phen$trans<- as.factor(phen$trans)
        w = (phen$tick_adult_M + phen$tick_adult_F)
phen = cbind(phen,w)        


        dim(cbind(phen[phen$year.y==2011,]$tick_nymph,phen[phen$year.y==2011,]$w))
     q = (phen[phen$year.y==2011,]$Temp)
        names(phen)
        
        
        
        dim(phen[phen$year.y==2011,]$Temp)
        

summary(M23B)                                           
plot((phen$tot)~(phen$Temp))
curve(((M23B$b[1]+1) + M23B$b[3]*x), from = 0, to = 10, add=T)                                           
curve(((M23B$b[1] + M23B$b[2]+) + (M23B$b[1] + M23B$b[6])*x), from = log(10), to = 30, add=T)                                    
load("25.11.climateworkspace.RData") 

nrow(phen)
nrow(phenuc)


plot(phen$temp, phen$tot)


dev.copy2eps(device=x11,
             file = "tull.eps",
             #paper="special",
             #width=10,
             #height=8,
             #horizontal=FALSE
)


save(list=ls(all=T), file="~/uio/phenology/01.12.climateworkspace.RData")