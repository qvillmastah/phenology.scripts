#cd uio/phenology
#R
rm(list = ls(all= TRUE))
#source("script.R", echo=T)
#install.packages("date")
library(date)

#################################################################################3333
########3333  Loading data

phen11 = read.table("dataflatt2011.csv", header=T, sep=",")
phen12 = read.table("dataflatt2012.csv", header=T, sep=",")
c = read.table("climate.csv", header=T, sep=",")


######  Transforing the phenolgy data
as.numeric(mdy.date(4, 26, 2011))-18628
as.numeric(mdy.date(11, 23, 2011))
attach(phen11)
date = as.numeric(mdy.date(month, day, year)-18628)
time.int = date*48 + hour*2
halvtime = min
head(min)
halvtime =as.numeric( ifelse(halvtime > 29, 1, 0))
head(halvtime) 
time.int = time.int+halvtime
detach(phen11)
phen11 = cbind(phen11, date, time.int)

attach(phen12)
date = as.numeric(mdy.date(month, day, year)-18627)
time.int = date*48 + as.numeric(hour)*2
halvtime = min
head(min)
halvtime =as.numeric( ifelse(halvtime > 29, 1, 0))
head(halvtime) 
time.int = time.int+halvtime
detach(phen12)
phen12 = cbind(phen12, date, time.int)


phenuc <- rbind(phen11, phen12)


######### Transforming climate data
attach(c)
date = as.numeric(mdy.date(month, day, year)-18627)
time.int = date*48 + hour*2
halvtime = min
head(min)
halvtime = as.numeric(ifelse(halvtime == 30, 1, 0) )
time.int = time.int+halvtime
SVP = 6.112*exp(17.67*Temp/(Temp+243.5))
VPD  = (100-Hum)*(SVP/100) # According to the formula by Gilbert 2010, SVP is calculated in hPa, while VPD uses kPa. But I get strange results. Seems like it should be divided by 100 instead of 10 anyway.
clim = data.frame(Temp, Hum, Dew, Place, date, SVP, VPD, time.int, year)
nrow(clim)
detach(c)
names(clim)
head(clim)
rm(VPD, SVP, time.int, halvtime)


clam = clim
clim =na.omit(clim)
clim$date = as.numeric(clim$date)
#clim$tot = as.numeric(clim$tot)
nrow(clim)
nrow(clam)

#############################################################################################3
####3   Checking out the climate data data, are any of the plots from the office?
############################################################################################333######333

clim11=subset(clim, clim$year==2011)
clim12=subset(clim, clim$year==2012)

is1_11 = subset(clim11, clim11$Place == "IsfjordenH1")
is2_11 = subset(clim11, clim11$Place == "IsfjordenH2")
tor1_11 = subset(clim11, clim11$Place == "TorjulvagenH1")
tor2_11 = subset(clim11, clim11$Place == "TorjulvagenH2")
is1_12 = subset(clim12, clim12$Place == "IsfjordenH1")
is2_12 = subset(clim12, clim12$Place == "IsfjordenH2")
tor1_12 = subset(clim12, clim12$Place == "TorjulvagenH1")
tor2_12 = subset(clim12, clim12$Place == "TorjulvagenH2")

#par(mfrow=c(3,3))
#plot(is1_11$Temp ~ is1_11$time.int, type = "n")#, xlim =c(0,1000))
#lines(is1_11$Temp ~ is1_11$time.int, col = "blue")

#plot(is2_11$Temp ~ is2_11$time.int, type = "n")#, xlim =c(15000,16000))
#lines(is2_11$Temp ~ is2_11$time.int, col = "blue")

#plot(is1_12$Temp ~ is1_12$time.int, type = "n")#, xlim =c(32000,33000))
#lines(is1_12$Temp ~ is1_12$time.int, col = "blue")

#plot(is2_12$Temp ~ is2_12$time.int, type = "n")#, xlim =c(32000,34000))
#lines(is2_12$Temp ~ is2_12$time.int, col = "blue")

#plot(tor1_11$Temp ~ tor1_11$time.int, type = "n")#, xlim =c(5000,8000))
#lines(tor1_11$Temp ~ tor1_11$time.int, col = "blue")

#plot(tor2_11$Temp ~ tor2_11$time.int, type = "n")#, xlim =c(6500,8000))
#lines(tor2_11$Temp ~ tor2_11$time.int, col = "blue")

#plot(tor1_12$Temp ~ tor1_12$time.int, type = "n")#, xlim =c(22000,25000))
#lines(tor1_12$Temp ~ tor1_12$time.int, col = "blue")

#plot(tor2_12$Temp ~ tor2_12$time.int, type = "n")#, xlim =c(32800,33200))
#lines(tor2_12$Temp ~ tor2_12$time.int, col = "blue")



is1_11 = is1_11[is1_11$time.int<15400&is1_11$time.int>7000,]
is2_11 = is2_11[is2_11$time.int<15400&is2_11$time.int>7000,]
is1_12 = is1_12[is1_12$time.int<33000&is1_12$time.int>24000,]
is2_12 = is2_12[is2_12$time.int<33000&is2_12$time.int>25000,]
tor1_11 = tor1_11[tor1_11$time.int<15400,]
tor2_11 = tor2_11[tor2_11$time.int<15400&tor2_11$time.int>7000,]
tor1_12 = tor1_12[tor1_12$time.int<32500&tor1_12$time.int>23500,]
tor2_12 = tor2_12[tor2_12$time.int<33000&tor2_12$time.int>24000,]


climate=rbind(is1_11,is2_11, is1_12,is2_12,tor1_11, tor2_11,tor1_12, tor2_12 )
#attach(climate)
#cor.test(Hum, Temp, "pearson", alternative = "two.sided")
#cor.test(Hum, Temp, "kendall", alternative = "two.sided")
#cor.test(Hum, Temp, "spearman", alternative = "two.sided")
#cor.test(VPD, Temp, "pearson", alternative = "two.sided")
#cor.test(VPD, Temp, "kendall", alternative = "two.sided")
#cor.test(VPD, Temp, "spearman", alternative = "two.sided")
#cor.test(VPD, Hum, "pearson", alternative = "two.sided")
#cor.test(VPD, Hum, "kendall", alternative = "two.sided")
#cor.test(VPD, Hum, "spearman", alternative = "two.sided")
#detach(climate)
#plot(climate$time.int, climate$Temp)
#plot(climate$time.int, climate$Hum)
#plot(climate$time.int, climate$VPD)




#attach(clim)
#a = aggregate(clim, by = list(date, as.character(Place)), FUN = mean, na.rm = T)
#head(a)
#detach(clim)
#nrow(a)


## Remember that H1 is low, H2 is high
#i.eb.evapfr -m


nrow(phen11)
nrow((phen11[phen11$Sted=="T"& phen11$alt.cat=="low",]))
nrow((phen11[phen11$Sted=="T"& phen11$alt.cat=="high",]))
nrow((phen11[phen11$Sted=="I"& phen11$alt.cat=="low",]))
nrow((phen11[phen11$Sted=="I"& phen11$alt.cat=="high",]))

isf.low11 <- merge(phen11[phen11$Sted=="I"& phen11$alt.cat=="low",],is1_11 ,by="time.int")
isf.high11 <- merge(phen11[phen11$Sted=="I"& phen11$alt.cat=="high",],is2_11 ,by="time.int")
tor.low11 <- merge(phen11[phen11$Sted=="T"& phen11$alt.cat=="low",],tor1_11 ,by="time.int")
tor.high11 <- merge(phen11[phen11$Sted=="T"& phen11$alt.cat=="high",],tor2_11 ,by="time.int")

nrow(isf.low11)
nrow(isf.high11)
nrow(tor.low11)
nrow(tor.high11)


nrow(phen12)
nrow((phen12[phen12$Sted=="T"& phen12$alt.cat=="low",]))
nrow((phen12[phen12$Sted=="T"& phen12$alt.cat=="high",]))
nrow((phen12[phen12$Sted=="I"& phen12$alt.cat=="low",]))
nrow((phen12[phen12$Sted=="I"& phen12$alt.cat=="high",]))

isf.low12 <- merge(phen12[phen12$Sted=="I"& phen12$alt.cat=="low",],is1_12 ,by="time.int")
isf.high12 <- merge(phen12[phen12$Sted=="I"& phen12$alt.cat=="high",],is2_12 ,by="time.int")
tor.low12 <- merge(phen12[phen12$Sted=="T"& phen12$alt.cat=="low",],tor1_12 ,by="time.int")
tor.high12 <- merge(phen12[phen12$Sted=="T"& phen12$alt.cat=="high",],tor2_12 ,by="time.int")

nrow(isf.low12)
nrow(isf.high12)
nrow(tor.low12)
nrow(tor.high12)

phen = rbind(isf.low11, isf.low12, isf.high11, isf.high12, tor.low11, tor.low12, tor.high11, tor.high12)
phen$year.x = as.factor(phen$year.x)



#db.out.ogr gps format=CSV dsn=/home/lars/gps_out.csv
#r.in.gdal input=/home/lars/MODIS/MODIS_250m_16_days_NDVI_Norden_without_h18v1_h19v1_2010_12_19_to_2011_01_03_EPSG_32633.tif output=kart2
#install.packages("AED", repos="http://www.highstat.com/book2.htm
#",type="source")


#r.sun --overwrite --verbose elevin=kart aspin=aspect2 slopein=slope day=150 insol_time=insol glob_rad=glob.rad #diff_rad=diff.rad refl_rad=refl.rad

library(splines)
library(glmmADMB)
library(MASS)
#library(VGAM)
## Enkel modell for å se om det er noen sammenheng mellom klimavariable og abundans

