#
source("data.R", echo=T)  # a separate script to load all the data and prepare
setwd("~/uio/phenology")
#load("~/uio/phenology/phenologyworkspace.RData")

phenuc$trans = as.factor(phenuc$trans)
head(phenuc)

AIC.list = vector(mode= "double", length=100)
i=1
## Staritng with only date
## 2011 - data
AIC.list[1]= AIC(mod1 <- glmmadmb(tot ~ ns(date, df = 1) + (1|uniktrans), data=phenuc[phenuc$year==2011,], family="nbinom"))
AIC.list[2]= AIC(mod2 <- glmmadmb(tot ~ ns(date, df=2) + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom"))
AIC.list[3]= AIC(mod3 <- glmmadmb(tot ~ ns(date, df=3) + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom"))
AIC.list[4]= AIC(mod4 <- glmmadmb(tot ~ ns(date, df=4) + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom"))
AIC.list[5]= AIC(mod5 <- glmmadmb(tot ~ ns(date, df=5) + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom"))

#2012 - data
AIC.list[7]= AIC(mod6 <- glmmadmb(tot ~ ns(date, df=1) + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom"))
AIC.list[8]= AIC(mod7 <- glmmadmb(tot ~ ns(date, df=2) + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom"))
AIC.list[9]= AIC(mod8 <- glmmadmb(tot ~ ns(date, df=3) + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom"))
AIC.list[10]= AIC(mod9 <- glmmadmb(tot ~ ns(date, df=4) + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom"))
AIC.list[11]= AIC(mod10 <- glmmadmb(tot ~ ns(date, df=5) + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom"))
       
## Including elevation. First for year == 2011
## First as additive effect
AIC.list[23]= AIC(mod23 <- glmmadmb(tot ~ ns(date, df=5) + alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
AIC.list[24]= AIC(mod24 <- glmmadmb(tot ~ ns(date, df=4) + alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
AIC.list[25]= AIC(mod25 <- glmmadmb(tot ~ ns(date, df=3) + alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
## Then as interaction
AIC.list[30]= AIC(mod30 <- glmmadmb(tot ~ ns(date, df=5) * alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
AIC.list[31]= AIC(mod31 <- glmmadmb(tot ~ ns(date, df=4) * alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
AIC.list[32]= AIC(mod32 <- glmmadmb(tot ~ ns(date, df=3) * alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2011,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))

## Including elevation. First for year == 2012
## First as additive effect
AIC.list[37]= AIC(mod37 <- glmmadmb(tot ~ ns(date, df=5) + alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
AIC.list[38]= AIC(mod38 <- glmmadmb(tot ~ ns(date, df=4) + alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
AIC.list[39]= AIC(mod39 <- glmmadmb(tot ~ ns(date, df=3) + alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE)))
## Then as interaction
AIC.list[44]= AIC(mod44 <- glmmadmb(tot ~ ns(date, df=5) * alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = TRUE, noinit = TRUE)))
AIC.list[45]= AIC(mod45 <- glmmadmb(tot ~ ns(date, df=4) * alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom", zeroInflation =T,admb.opts = admbControl(shess = FALSE, noinit = FALSE))) 
AIC.list[46]= AIC(mod46 <- glmmadmb(tot ~ ns(date, df=3) * alt.cat + (1|Sted/trans), data=phenuc[phenuc$year==2012,], family="nbinom", admb.opts = admbControl(shess = FALSE, noinit = FALSE)))


AIC.list

## Plotting the models. Here you can just change the model in object "model best - and see what you get"
   #              
par(mfrow=c(2,1))                 
       
model.best = mod30  ## The best AIC that does not add arbritary trends to the predictions - 2011
newdat <- expand.grid(date=c(117:326), alt.cat=c("high","low"))
mm <- model.matrix(delete.response(terms(model.best)),newdat)
newdat$tot <- exp(mm %*% fixef(model.best))
predvar <- diag(mm %*% vcov(model.best) %*% t(mm))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+model.best$alpha^2)
newdat$tot=exp(newdat$tot)
newdat$low=exp(newdat$tot-newdat$SE2)
newdat$high=exp(newdat$tot+newdat$SE2)



plot(newdat$date, log(newdat$tot), type="n",  ylim=c(0,4), main ="2011", ylab = "tick counts", xlab="", xaxt = "n", las = 1, bty="n")
points(phenuc[phenuc$year==2011&phenuc$alt.cat=="low",]$date, log(phenuc[phenuc$year==2011&phenuc$alt.cat=="low",]$tot), pch=25, bg="gray70", col="gray70",cex=1)
points(phenuc[phenuc$year==2011&phenuc$alt.cat=="high",]$date, log(phenuc[phenuc$year==2011&phenuc$alt.cat=="high",]$tot), pch=2, col="black",cex=1)
attt = c(90,120,151, 181, 212, 243, 273, 304, 334)
attt = attt
#labs = as.date(attt+ 18628)
labs=c("Apr","May","Jun","Jul","Aug","Sep","Oct", "Nov", "Dec")
axis(1, at = attt, labels = labs, col="black")
lines(newdat$date[newdat$alt.cat == "low"&newdat$date>117&newdat$date<326], log(newdat$tot[newdat$alt.cat == "low"&newdat$date>117&newdat$date<326]), lty =2, lwd=2)
lines(newdat$date[newdat$alt.cat == "high"&newdat$date>156&newdat$date<299], log(newdat$tot[newdat$alt.cat == "high"&newdat$date>156&newdat$date<299]), lwd=2)

#                 legend(attt[8], 4.2, c("Low", "High", "Low", "High"),
#        text.col = "black",
#        lty = c(NA, NA, 2,1), 
#        pch = c(25, 2, NA,NA),
#	col = c("gray70", "black", "black","black"),
#        pt.bg = c("gray70",NA, NA, NA),
#       bty="n", bg ="white"
#)
                 
                 
                 
model.best = mod44  ## The best AIC that does not add arbritary trends to the predictions - 2012
newdat <- expand.grid(date=c(482:662), alt.cat=c("high","low"))
mm <- model.matrix(delete.response(terms(model.best)),newdat)
newdat$tot <- exp(mm %*% fixef(model.best))
predvar <- diag(mm %*% vcov(model.best) %*% t(mm))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+model.best$alpha^2)
newdat$tot=exp(newdat$tot)
newdat$low=exp(newdat$tot-newdat$SE2)
newdat$high=exp(newdat$tot+newdat$SE2)
                 
plot(newdat$date, log(newdat$tot), type="n",  ylim=c(0,4), main ="2012", ylab = "tick counts", xlab="", xaxt = "n", las = 1, bty="n")
points(phenuc[phenuc$year==2012&phenuc$alt.cat=="low",]$date, log(phenuc[phenuc$year==2012&phenuc$alt.cat=="low",]$tot), pch=25, bg="gray70", col="gray70",cex=1)
points(phenuc[phenuc$year==2012&phenuc$alt.cat=="high",]$date, log(phenuc[phenuc$year==2012&phenuc$alt.cat=="high",]$tot), pch=2, col="black",cex=1)
attt = c(90,120,151, 181, 212, 243, 273, 304)
attt = attt+365
labs=c("Apr","May","Jun","Jul","Aug","Sep","Oct", "Nov")
axis(1, at = attt, labels = labs, col="black")
lines(newdat$date[newdat$alt.cat == "low"&newdat$date>482&newdat$date<662], log(newdat$tot[newdat$alt.cat == "low"&newdat$date>482&newdat$date<662]), lty =2, lwd=2)
lines(newdat$date[newdat$alt.cat == "high"&newdat$date>509&newdat$date<649], log(newdat$tot[newdat$alt.cat == "high"&newdat$date>509&newdat$date<649]), lwd=2)

                 legend(attt[7], 4, c("Low", "High", "Low", "High"),
        text.col = "black",
        lty = c(NA, NA, 2,1), 
        pch = c(25, 2, NA,NA),
	col = c("gray70", "black", "black", "black"),
        pt.bg = c("gray70",NA, NA, NA),
        bty="n",
        bg = "white"                        
)
       
#       
dev.copy2eps(device=x11,
file = "~/uio/phenology/figures/Figure1.eps",
paper="special",
    #   width=10,
    #   height=10,
      # horizontal=FALSE
		)
dev.copy2pdf(device=x11,
file = "~/uio/phenology/figures/Figure1.pdf",
paper="special",
     #  width=10,
    #   height=10,
     #  horizontal=FALSE
  	)



## Onset
#Spring 2011, low elevation
min(phenuc[phenuc$year==2011,]$date)
phenuc[phenuc$date==117,]
#Spring 2012, low elevation
min(phenuc[phenuc$year==2012,]$date)
phenuc[phenuc$date==482,]
#Spring 2011, high elevation
min(phenuc[phenuc$year==2011&phenuc$alt.cat =="high",]$date)
phenuc[phenuc$date==156,]
#Spring 2012, high elevation
min(phenuc[phenuc$year==2012&phenuc$alt.cat =="high",]$date)
phenuc[phenuc$date==509,]


## End of questing season
#Spring 2011, low elevation
max(phenuc[phenuc$year==2011,]$date)
phenuc[phenuc$date==326,]
#Spring 2012, low elevation
max(phenuc[phenuc$year==2012,]$date)
phenuc[phenuc$date==662,]
#Spring 2011, high elevation
max(phenuc[phenuc$year==2011&phenuc$alt.cat =="high",]$date)
phenuc[phenuc$date==299,]
#Spring 2012, high elevation
max(phenuc[phenuc$year==2012&phenuc$alt.cat =="high",]$date)
phenuc[phenuc$date==649,]

dev.off()
           
                 
coplot(tot~ date| alt.cat + as.factor(year), data = phenuc)


mod2011 <- mod31
mod2012 <- mod44
save(mod2011, file="mod2011.RData")
save(mod2012, file="mod2012.RData")
save(list=ls(all=T), file="~/uio/phenology/phenologyworkspace.RData")