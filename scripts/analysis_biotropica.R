# Load libraries ----
library(tidyverse)
library(lme4)
library(reshape)
library(nlme) 
library(MuMIn)
library(AICcmodavg)

# Load data ----

source("./scripts/dataPrep.R")

# Litter quality ----

t.test(CN ~ Source, data=chns)
cnmod <- lmer(CN ~ Source + (1|Batch), data=chns)
anova(cnmod, update(cnmod, ~.-Source))

# DECOMPOSITION ----

fullmodel <- lme(LdryMass ~ LinitialMass + 
                   numDaysField +
		numDaysField:Site_Cond + numDaysField:Source_Cond +
		numDaysField:Site_Cond:Source_Cond,
		random=list(~1|Site, ~1|Batch), 
		weights= varExp(form =~ as.vector(numDaysField) ),	
		data=decomp, subset=complete.cases(decomp), 
		na.action=na.fail, method="ML")
plot(fullmodel)
(MMItab <- dredge(fullmodel, trace=T, subset=LinitialMass))

# refit with REML
bestmod <- lme(LdryMass ~ LinitialMass + 
                 numDaysField +
                 numDaysField:Site_Cond + numDaysField:Source_Cond +
                 numDaysField:Site_Cond:Source_Cond,
               random = list(~1 | Site, ~1 | Batch),
               weights = varExp(form = ~as.vector(numDaysField)),
               data = decomp, subset = complete.cases(decomp),
               na.action = na.fail, method = "REML")
summary(bestmod)

## MODEL PREDICTIONS ----

# new day variable
newDays<- 0:300
newdat.wet <- data.frame(Site_Cond="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)))
newdat.dry <- data.frame(Site_Cond="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)))

# predictions
dry.pred <- predictSE(bestmod, newdata=newdat.dry, level=0, se.fit=T)
wet.pred <- predictSE(bestmod, level=0, newdata=newdat.wet, se.fit=T)

bmc <- summary(bestmod)$coef$fixed
dry.pred <- with(newdat.dry, bmc[1] + LinitialMass*bmc[2] + numDaysField*bmc[3] 
	+ ifelse(Site_Cond=="wet",1,0)*numDaysField*bmc[4])
wet.pred <- with(newdat.wet, bmc[1] + LinitialMass*bmc[2] + numDaysField*bmc[3] 
	+ ifelse(Site_Cond=="wet",1,0)*numDaysField*bmc[4])

cols <- c("#70A37F","#41658A")
cols <- c("black", "grey30", "white")
colt <- adjustcolor(cols, alpha.f=0.2)

summary(decomp)
head(decomp)
table(decomp$Batch, decomp$Site)
table(decomp$fill.date, decomp$Batch)

## BOOTSTRAPPING ---- 
### SPLIT ESTIMATES BY ORIGIN ----

dryplots <- c("213", "408", "107", "405", "1")
wetplots <- c("6", "9", "10", "3", "4")

diff25 <- list()

for(i in 1:5){
  for(j in 1:5){
    diff25[[5*(i-1)+j]] <-
      with(subset(decomp, Site==dryplots[i]), tapply(propReduced, list(batch, months, Source_Cond), mean, na.rm=T)) -
      with(subset(decomp, Site==wetplots[j]), tapply(propReduced, list(batch, months, Source_Cond), mean, na.rm=T))
  }  }

diff25[[25]][,,"dry"]

drysource <- matrix(rep(NA,15), nrow=5)
wetsource <- matrix(rep(NA,15), nrow=5)
colnames(drysource) <- c("lower", "median", "upper")
colnames(wetsource) <- c("lower", "median", "upper")
m <- c(0,1,2,4,8)
rownames(drysource) <- m
rownames(wetsource) <- m

for(i in 1:5){
  
  drysource[i,] <- 
    quantile(
      sample(
        na.omit(unlist(lapply(diff25, "[" , ,i,"dry")))
        , 9999, replace=T)
      , c(0.025, 0.5, 0.975))
  
  wetsource[i,] <- 
    quantile(
      sample(
        na.omit(unlist(lapply(diff25, "[" , ,i,"wet")))
        , 9999, replace=T)
      , c(0.025, 0.5, 0.975))
}

par(mar=c(6.5,6,4.5,2), mgp=c(4,1,0))
plot(1~1, xlim=c(-0.5,8.5), ylim=range(rbind(wetsource,drysource)),
     xlab="Number of months in the field", ylab="Difference in decomposed litter proportion",
     cex.lab=2, cex.axis=1.5, las=1)
rect(-1,-1,9,0, border=F, col="grey95")
abline(h=0, lty=2, lwd=2)
box()
for(i in 1:5){
  arrows(c(m-0.15)[i], wetsource[i,"lower"], c(m-0.15)[i], wetsource[i,"upper"],
         length=0, lwd=3)}
for(i in 1:5){
  arrows(c(m+0.15)[i], drysource[i,"lower"], c(m+0.15)[i], drysource[i,"upper"],
         length=0, lwd=3, col="grey60")}
points(wetsource[,"median"] ~ c(m-0.15), cex=3, pch=21, lwd=3, bg="white")
points(drysource[,"median"] ~ c(m+0.15), pch=21, cex=3, lwd=2, bg="grey60", col="white")
text(c(5.5,5), c(0.05,-0.025), cex=1.4, 
     labels=c("Litter in non-swamp plots\n decomposed more", "Litter in swamp plots decomposed more"))

### SPLIT ESTIMATES BY MICROHABITAT/PLOT OF DECOMPOSITION ----

dsite <- with(subset(decomp, Source_Cond=="dry"), tapply(propReduced, list(batch, months, Site), mean, na.rm=T)) - 
  with(subset(decomp, Source_Cond=="wet"), tapply(propReduced, list(batch, months, Site), mean, na.rm=T))
dmwd <- rbind(dsite[,,"3"],dsite[,,"4"],dsite[,,"6"],dsite[,,"9"],dsite[,,"10"])
dmww <- rbind(dsite[,,"1"],dsite[,,"213"],dsite[,,"405"],dsite[,,"408"],dsite[,,"107"])

drysource <- matrix(rep(NA,15), nrow=5)
wetsource <- matrix(rep(NA,15), nrow=5)
colnames(drysource) <- c("lower", "median", "upper")
colnames(wetsource) <- c("lower", "median", "upper")
m <- c(0,1,2,4,8)
rownames(drysource) <- m
rownames(wetsource) <- m

for(i in 1:5){
  
  drysource[i,] <- 
    quantile(
      sample(
        na.omit(dmwd[,i])
        , 9999, replace=T)
      , c(0.025, 0.5, 0.975))
  
  wetsource[i,] <- 
    quantile(
      sample(
        na.omit(dmww[,i])
        , 9999, replace=T)
      , c(0.025, 0.5, 0.975))
}

par(mar=c(6.5,6,4.5,2), mgp=c(4,1,0))
plot(1~1, xlim=c(-0.5,8.5), ylim=range(rbind(wetsource,drysource)),
     xlab="Number of months in the field", ylab="Difference in decomposed litter proportion",
     cex.lab=2, cex.axis=1.5, las=1)
rect(-1,-1,9,0, border=F, col="grey95")
abline(h=0, lty=2, lwd=2)
box()
for(i in 1:5){
  arrows(c(m-0.15)[i], wetsource[i,"lower"], c(m-0.15)[i], wetsource[i,"upper"],
         length=0, lwd=3)}
for(i in 1:5){
  arrows(c(m+0.15)[i], drysource[i,"lower"], c(m+0.15)[i], drysource[i,"upper"],
         length=0, lwd=3, col="grey60")}
points(wetsource[,"median"] ~ c(m-0.15), cex=3, pch=21, lwd=3, bg="white")
points(drysource[,"median"] ~ c(m+0.15), pch=21, cex=3, lwd=2, bg="grey60", col="white")
text(c(5,5), c(0.3,-0.2), cex=1.5, 
     labels=c("Non-swamp origin litter\n decomposed more", "Swamp origin litter\ndecomposed more"))

## PLOT ----

### SUM OF WEIGHTS ----

sow <- matrix(rep(NA,6))
rownames(sow) <- names(MMItab)[1:6]
for(i in 1:6){
	sow[i,1] <- sum(MMItab$weight[which(!is.na(MMItab[,i]))])}
(sow <- sow[order(sow[,1]),])
names(sow) <- c("Duration x Plot x Origin", "Duration x Origin",
	"Handling loss", "Initial mass", "Duration", "Duration x Plot")

#pdf("Fig 2 Oeco.pdf", height=15, width=6.5)
#jpeg("Fig 2 Oeco.jpg", height=15, width=7, res=300, units="in")
#setEPS(); postscript("Fig 2 Oeco.eps", width=6.5, height=15)
layout(matrix(1:3, ncol=1), heights=c(1,2,2))
	# a #
par(mar=c(5.5,20,4,4))
barplot(sow, beside=T, horiz=T, las=1, xlab="Sum of variable weights",
	cex.lab=1.5, cex.axis=1.5, cex.names=1.5, tck=0.02)
mtext(side=3, line=-2.5, text=" a) Sum of variable weights", cex=1.5, adj=0.2, outer=T)

	# b #
par(mar=c(6.5,6,2.5,2), mgp=c(4,1,0))
plot(dryMass ~ jitter(numDaysField,10), data = decomp, type="n",
	cex.lab=2, cex.axis=1.5, xlim=range(newDays), yaxt="n", tck=0.01,
	ylab="Litter dry mass (g)",  xlab="Number of days in field")
#points(dryMass ~ jitter(numDaysField,20), data = decomp, 
#	subset=Site_Cond=="dry", pch=21, col="grey", lwd=3, bg="white", cex = 3)
#points(dryMass ~ jitter(numDaysField,20), data = decomp, 
#	subset=Site_Cond=="wet", pch=16, col=colt[1], cex = 3)
points(dryMass ~ jitter(numDaysField,20), data = decomp,
	subset=Site_Cond=="dry", pch=21, col="grey60", lwd=1, bg="white", cex = 2)
points(dryMass ~ jitter(numDaysField,20), data = decomp, 
	subset=Site_Cond=="wet", pch=21, col="white", bg="black", lwd=1, cex = 2)
axis(2, at=0:6, labels=0:6, cex.axis=1.5, las=1, tck=0.01)

#lines(exp(dry.pred$fit) ~ newDays, col="white", lwd=12)
#lines(exp(wet.pred$fit) ~ newDays, col="white", lwd=12)
lines(exp(wet.pred) ~ newDays, col=cols[1], lwd=8)
lines(exp(dry.pred) ~ newDays, col="grey60", lwd=8, lty=2)

legend("topright", bty="n", cex=2,
	lty=c(1,2), lwd=4, col=c("black","grey50"),
	legend=c("Swamp","Non-swamp"), title="Plot type")
legend("topright", bty="n", cex=2,
	pt.cex=2, col=c("white","grey60"), pch=21, pt.bg=c("black","white"),
	legend=c("                        ",""), title="               ")
mtext(side=3, line=1, text=" b) Best model prediction", cex=1.5, adj=0)

	# c #
par(mar=c(6.5,6,4.5,2))
plot(c(1)~1, xlim=c(-0.5,8.5), ylim=range(rbind(wetsource,drysource)),
	xlab="Number of months in the field", ylab="Difference in decomposed litter proportion",
	cex.lab=2, cex.axis=1.5, las=1, tck=0.01)
mtext(side=3, line=1, text=" c) Effects of origin and plot type on decomposition", cex=1.4, adj=0)
rect(-1,-1,9,0, border=F, col="grey95")
abline(h=0, lty=2, lwd=2)
axis(side=1, at=seq(0,8,2), tck=0.01, labels=NA)
axis(side=2, at=seq(-0.2,0.3,0.1), tck=0.01, labels=NA)
box()
for(i in 1:5){
  arrows(c(m-0.15)[i], wetsource[i,"lower"], c(m-0.15)[i], wetsource[i,"upper"],
    length=0, lwd=3)}
for(i in 1:5){
  arrows(c(m+0.15)[i], drysource[i,"lower"], c(m+0.15)[i], drysource[i,"upper"],
    length=0, lwd=3, col="grey50")}
points(wetsource[,"median"] ~ c(m-0.15), cex=3, pch=21, lwd=3, bg="black", col="white")
points(drysource[,"median"] ~ c(m+0.15), pch=21, cex=3, lwd=3, bg="white", col="grey50")
text(c(4,4), c(0.3,-0.2), cex=1.5, 
	labels=c("Non-swamp origin litter\n decomposed more", "Swamp origin litter\ndecomposed more"))
legend("topright", title="Plot type", legend=c("Non-swamp","Swamp"), 
	pch=21, cex=1.5, bty="n", pt.cex=3, col=c("white","grey50"), pt.bg=c("black","white"), pt.lwd=2)
	
dev.off()

### SOURCE VS SITE BARPLOTS ----

eight <- droplevels(decomp[which(decomp$numDaysField>1.5),])
means <- with(na.omit(eight), tapply(dryMass, list(Source_Cond, Site_Cond), mean))
ses <- with(eight, tapply(dryMass, list(Source_Cond, Site_Cond), se))
low <- c(means-ses)
upp <- c(means+ses)

bp <- barplot(means, beside=T, col=cols, ylim=c(0,4), 
	xlab="Plot type", ylab="Litter dry mass after 8 months (g)",
	cex.lab=2, cex.names=2, cex.axis=1.5, las=1)
for(i in 1:4){
	arrows(bp[i], low[i], bp[i], upp[i], length=0.25, angle=90, code=3, lwd=1)}
legend(0.8, 3.7, bty="n", fill=cols, title="Litter origin", 
	legend=c("Non-swamp tree species", "Swamp tree species"), cex=2)
mtext(side=3, line=-0.5, text=" c) Remaining litter after 8 months", cex=1.5, adj=0)

dev.off()
with(decomp, tapply(propRemain, list(months, Site_Cond), mean, na.rm=T))

# Table 1 ----

### Get all top models plus best model containing HFA:Duration
topmodels <- get.models(MMItab, subset=T)
t2<-length(topmodels)
top2<-vector("list",length=t2)

mm <- MMItab

top.coefs <- data.frame(
	Handling.loss=rep(NA,t2),
	Initial.mass=rep(NA,t2),
	Duration=rep(NA,t2),
	DurationxPlot=rep(NA,t2),
	DurationxOrigin=rep(NA,t2),
	DurationxPlotxOrigin=rep(NA,t2))
var.names = c("(Intercept)", "LinitialMass", "numDaysField", 
	"numDaysField:Site_Condwet", "numDaysField:Source_Condwet",
	"numDaysField:Site_Condwet:Source_Condwet")
names(top.coefs);var.names

for(j in 1:length(top.coefs)){
	for(i in 1:t2){	
		tryCatch({
		  top.coefs[i,j] <- paste0(
		    round(
			summary(topmodels[[i]])$tTable[var.names[j],1],
		    4),
		  " (?",
		    round(
		      1.96*summary(topmodels[[i]])$tTable[var.names[j],2],
		    4),
		  ")")
		},error=function(e){})		
	}
}
top.coefs
a <- data.frame(round(mm[,10:11], 3))
Rsquare<-data.frame(R2m=rep(0,t2),R2c=rep(0,t2))
for(i in 1:length(topmodels)){
	Rsquare[i,]<-round(r.squaredLR(topmodels[[i]]),3)}
(table <- data.frame(top.coefs, df=mm$df, a, Rsquare))
#write.csv(table, "Model selection table v2.csv")      

# BEST MODEL COEF PLOT ----

bm.coef <- data.frame(summary(bestmod)$tTable[-1,1:2])
bm.coef$upp <- bm.coef$Value + 1.96*bm.coef$Std.Error
bm.coef$low <- bm.coef$Value - 1.96*bm.coef$Std.Error

par(mar=c(5,15,2,2))
plot(1:4 ~ Value, xlim=c(min(low), max(upp)), data=bm.coef, pch=16, cex=4,
	cex.lab=2, cex.axis=1.5, yaxt="n", ylab="", xlab="Effect size",
	ylim=c(0.5,4.5))
for(i in 1:4){
	arrows(bm.coef$upp[i], i, bm.coef$low[i], i, length=0, lwd=5)}
abline(v=0, lty=2, lwd=2)
axis(side=2, at=1:4, cex.axis=1.5, las=1,
	lab=c("Initial mass", "Duration", "Site condition = wet", "Duration x\nSite condition = wet"))

### MODEL AVERAGING ----
MMIavg <- model.avg(subset(MMItab, delta<2))
summary(MMIavg)
am.coef <- data.frame(summary(MMIavg)$coefmat.full[-1,c(1,3)])
am.coef$upp <- am.coef$Estimate + 1.96*am.coef[,2]
am.coef$low <- am.coef$Estimate - 1.96*am.coef[,2]
am.coef
am.coef <- am.coef[c(6,4,5,3,2,1),]

col.dot <- ifelse(am.coef$upp < 0 | am.coef$low > 0, "black", "grey")

par(mar=c(5,15,2,2))
plot(1:6 ~ Estimate, xlim=c(min(low), max(upp)), data=am.coef, ylim=c(0.5,6.5),
	cex.lab=2, cex.axis=1.5, yaxt="n", ylab="", xlab="Effect size")
abline(v=0, lty=2, lwd=2)
points(1:6 ~ Estimate, xlim=c(min(low), max(upp)), data=am.coef, pch=16, cex=4, col=col.dot)
for(i in 1:6){
	arrows(am.coef$upp[i], i, am.coef$low[i], i, length=0, lwd=5, col=col.dot[i])}
axis(side=2, at=1:6, cex.axis=1.5, las=1,
	lab=c("Home field advantage", "Duration x\nSite condition = wet", "Source condition = wet", 
	"Site condition = wet", "Duration", "Initial mass"))

pairs.cor <- function (x,y,smooth=TRUE, digits=2,  ...)
{
  panel.cor <- function(x, y, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r.obj = cor.test(x, y,use="pairwise",...)
    r = as.numeric(r.obj$estimate)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(txt)
    cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex=cex*abs(r))
  }
panel.hist <- function(x)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="grey")
  }
pairs(x,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.smooth, ...)
}
