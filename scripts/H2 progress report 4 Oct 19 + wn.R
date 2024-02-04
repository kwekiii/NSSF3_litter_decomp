
setwd("F:/National University of Singapore/Chong Kwek Yan - NSSF3/CRSF/Past projects - DO NOT SHARE/R_JJ")
decomp<-read.csv("decomp_30.9.csv", header = T)
# Load libraries
library(lme4)
library(reshape)
library(nlme) 
library(MuMIn)


#preparing data
decomp$Site<-as.factor(decomp$Site)
decomp$TimeFrame<-as.factor(decomp$TimeFrame)
numDaysField.unscaled<-decomp$numDaysField
decomp$numDaysField<-scale(decomp$numDaysField)
decomp$LdryMass<-log(decomp$dryMass)
decomp$LinitialMass<-log(decomp$initialMass)
summary(decomp)

##keep batach as ran eff##

fullmodel.batchrf <- lme(LdryMass ~ LinitialMass + numDaysField + Site_Cond + Source_Cond +
                         Source_Cond:numDaysField + Site_Cond:numDaysField + Source_Cond:Site_Cond +
                         Source_Cond:Site_Cond:numDaysField, 
                       random=list(~1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	
                       data=decomp, subset=complete.cases(decomp), na.action=na.fail, method="ML")
MMItab.batchrf<-dredge(fullmodel.batchrf, subset = !(Source_Cond))

#write.csv(MMItab.batchrf, "dredgeresult.csv", quote=FALSE)
bestmods.batchrf <- get.models(MMItab.batchrf, subset = delta<2)

# new day variable
newDays.unscaled <- 0:300
estDays <-  (newDays.unscaled - mean(numDaysField.unscaled)) / sd(numDaysField.unscaled)
xlabels <- seq(0,300,50)
xticks <- estDays[seq(0,300,50)+1]

newdat.wet <- data.frame(Site_Cond="wet",numDaysField=estDays,LinitialMass=mean(log(decomp$initialMass)))
newdat.dry <- data.frame(Site_Cond="dry",numDaysField=estDays,LinitialMass=mean(log(decomp$initialMass)))

# predictions
dry.pred <- predict(bestmods.batchrf[[1]], level=0, newdata=newdat.dry)
wet.pred <- predict(bestmods.batchrf[[1]], level=0, newdata=newdat.wet) 
summary(bestmods.batchrf[[1]])

# plot
#jpeg("Litter decomp.jpg", units="in", height=9, width=12, res=300)
par(mar=c(5,5,2,2))
plot(dryMass ~ jitter(numDaysField), data = decomp, 
	col = ifelse(decomp$Site_Cond=="wet", "#41658A80", "#70A37F80"), pch=16,
	cex = 2, cex.lab=1.5, cex.axis=1.3, xaxt="n", xlim=range(estDays),
	ylab="Litter dry mass (g)",  xlab="Number of days in field", yaxt="n")
p <- pretty(par("usr")[3:4])
l <- formatC(p,format="f", digits=2)
axis(2, at=p, labels=l, cex.axis=1.2)
axis(1, at=xticks, labels=xlabels, cex.axis=1.2)

legend("topright", bty="n", cex=1.2,
	lty=c(1,2),lwd=4,
	legend=c("Wet","Dry"), title="Plot condition",
       col=c("#41658AB3","#70A37FB3"), pch=16, pt.cex=2)

lines(exp(dry.pred) ~ estDays, col="#70A37F", lwd=5, lty=2)
lines(exp(wet.pred) ~ estDays, col="#41658A", lwd=5)
dev.off()


        
      