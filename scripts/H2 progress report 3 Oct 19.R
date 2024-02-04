
setwd("D:/National University of Singapore/Chong Kwek Yan - FYP 1819 Lian Jun Jie/R")
decomp<-read.csv("decomp_30.9.csv", header = T)
# Load libraries
library(lme4)
library(AICcmodavg)
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

initial.model<-lmer(log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
                       Ratio:SiteC:numDaysField +(1|Site)+ (1|Batch),	data=decomp,REML=T)
summary(initial.model)
plot(initial.model)

initial.model2<-lmer(log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
                      Ratio:SiteC:numDaysField +(1|Site),	data=decomp,REML=T)
plot(initial.model2)

anova(initial.model,initial.model2)
##keep batach as ran eff##

var.model1<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
                  Ratio:SiteC:numDaysField  , random=list(~ 1|Site, ~1|Batch), weights= varPower(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude)

var.model2<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
                  Ratio:SiteC:numDaysField          , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude)


var.model3<-lme(fixed=LdryMass ~ LinitialMass + numDaysField + Ratio + SiteC 
           + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField + Ratio:SiteC:numDaysField, 
           random=list(~ 1|Site,~1|Batch), weights= varConstPower(form =~ as.vector(numDaysField) ), data=decomp, na.action=na.exclude, method="ML")

AIC(var.model1,var.model2,var.model3)
#var model2 lowest AIC

rmodel1<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
                  Ratio:SiteC:numDaysField          , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude)

rmodel2<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
               Ratio:SiteC:numDaysField          , random=list(~ 1|Site), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude)

anova(rmodel1,rmodel2)
#with batch better
#but, batch is directly associated with Ratio...

model<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
               Ratio:SiteC:numDaysField          , random=list(~ 1|Site), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
plot(model)
summary(model)

model1<-update(model, ~.-Ratio:SiteC:numDaysField)
anova(model, update(model, ~.-Ratio:SiteC:numDaysField))

summary(model1)

model2a<-update(model1, ~.-Ratio:numDaysField)
model2b<-update(model1, ~.-Ratio:SiteC)

anova(model1,model2a)
anova(model1,model2b)

model3<-update(model2a, ~.-Ratio:SiteC)
anova(model2a,model3)
summary(model3)

model4<-update(model3,~.-Ratio)
anova(model3,model4)

model<-lme(fixed=log(dryMass) ~ LinitialMass + numDaysField + Ratio + SiteC + SiteC:Ratio + SiteC:numDaysField + Ratio:numDaysField+
             Ratio:SiteC:numDaysField          , random=list(~ 1|Site), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")


new.wet.data.RatioH<-data.frame(SiteC="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.dry.data.RatioH<-data.frame(SiteC="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.wet.data.RatioL<-data.frame(SiteC="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.dry.data.RatioL<-data.frame(SiteC="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))


plot(dryMass~numDaysField.unscaled,data=decomp, type="n")

points(dryMass~numDaysField.unscaled, data=decomp, subset = c(SiteC=="dry"), col="pink", pch=2, cex=Ratio/50)
points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="wet"), col="lightblue", pch=1, cex=Ratio/50)

lines(exp(dry.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=1)

lines(exp(dry.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=2)

lines(exp(wet.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=1)

lines(exp(wet.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=2)


  anova(model, update(model, ~.-Ratio:SiteC:numDaysField))
  
  
fullmodel.batchrf<-lme(LdryMass~LinitialMass+numDaysField+SiteC+Source+Ratio+
                         Source:numDaysField+Ratio:numDaysField+SiteC:numDaysField+
                         Source:SiteC:numDaysField+Ratio:SiteC:numDaysField, 
                       random=list(~1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	
                       data=decomp, subset=complete.cases(decomp), na.action=na.fail, method="ML")

fullmodel.batchrfb<-lme(LdryMass~LinitialMass+numDaysField+SiteC+Source+Ratio+
                         Source:numDaysField+Ratio:numDaysField+SiteC:numDaysField+
                         Source:SiteC:numDaysField+Ratio:SiteC:numDaysField, 
                       random=list(~1|Site, ~1|Batch), weights= varConstPower(form =~ as.vector(numDaysField) ),	
                       data=decomp, subset=complete.cases(decomp), na.action=na.fail, method="ML")  

#use varExp

  
fullmodel<-lme(LdryMass~
                 LinitialMass+numDaysField+SiteC+Source+Ratio+
                 Source:numDaysField+Ratio:numDaysField+SiteC:numDaysField+
                 Source:SiteC:numDaysField+Ratio:SiteC:numDaysField, 
               random=list(~1|Site), weights= varExp(form =~ as.vector(numDaysField) ),	
               data=decomp, subset=complete.cases(decomp), na.action=na.fail, method="ML")

MMItab.batchrf<-dredge(fullmodel.batchrf, subset = !(Source & Ratio))
MMItab<-dredge(fullmodel, subset = !(Source & Ratio))

write.csv(MMItab.batchrf, "dredgeresult.csv", quote=FALSE)
write.csv(MMItab, "MMItab.csv", quote=FALSE)

bestmods.batchrf<-get.models(MMItab.batchrf, subset = delta<2)

library(AICcmodavg)

newDays <- seq(min(decomp$numDaysField,na.rm=TRUE),max(decomp$numDaysField,na.rm=TRUE), length.out = 1000)
newDays.unscaled<-(newDays*sd(numDaysField.unscaled)+mean(numDaysField.unscaled))

new.dw.data.RatioH<-data.frame(SiteC="dry",Source="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.wd.data.RatioH<-data.frame(SiteC="wet",Source="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.ww.data.RatioH<-data.frame(SiteC="wet",Source="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.dd.data.RatioH<-data.frame(SiteC="dry",Source="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.ww.data.RatioL<-data.frame(SiteC="wet",Source="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.dw.data.RatioL<-data.frame(SiteC="dry",Source="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.wd.data.RatioL<-data.frame(SiteC="wet",Source="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.dd.data.RatioL<-data.frame(SiteC="dry",Source="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))

dd.pred.RatioH<-modavgPred(bestmods.batchrf, newdata=new.dd.data.RatioH)
ww.pred.RatioH<-modavgPred(bestmods.batchrf, newdata=new.ww.data.RatioH)
dw.pred.RatioH<-modavgPred(bestmods.batchrf, newdata=new.dw.data.RatioH)
wd.pred.RatioH<-modavgPred(bestmods.batchrf, newdata=new.wd.data.RatioH)
dd.pred.RatioL<-modavgPred(bestmods.batchrf, newdata=new.dd.data.RatioL)
ww.pred.RatioL<-modavgPred(bestmods.batchrf, newdata=new.ww.data.RatioL)
dw.pred.RatioL<-modavgPred(bestmods.batchrf, newdata=new.dw.data.RatioL)
wd.pred.RatioL<-modavgPred(bestmods.batchrf, newdata=new.wd.data.RatioL)

lines(exp(dd.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2)
lines(exp(dw.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=2)
lines(exp(ww.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2)
lines(exp(wd.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=2)

lines(exp(dd.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=1)
lines(exp(dw.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=1, lty=2)
lines(exp(ww.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=1)
lines(exp(wd.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=1, lty=2)

decomp.col<-with(decomp, ifelse(SiteC=="dry",
                                ifelse(Source=="dry", "DarkGreen", "Red"),
                                ifelse(Source=="dry", "Sienna", "DodgerBlue")
))

pch=with(decomp,ifelse(SiteC=='dry',ifelse(Source=="dry",2,17), ifelse(Source=="dry",1,16)))

par(mfrow=c(1,2))
plot(dryMass ~ numDaysField.unscaled, data = decomp, col = decomp.col, cex = 2, font.axis = 
       1, font.lab = 2, pch=pch, main="Title", ylab="Dry Mass",  xlab="Number of Days in Field", xlim=c(0,160), type="n")
lines(exp(dd.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2)
lines(exp(dw.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=2)
lines(exp(ww.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2)
lines(exp(wd.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=2)

plot(dryMass ~ numDaysField.unscaled, data = decomp, col = decomp.col, cex = 2, font.axis = 
       1, font.lab = 2, pch=pch, main="Title", ylab="Dry Mass",  xlab="Number of Days in Field", xlim=c(0,160), type="n")
lines(exp(dd.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=1)
lines(exp(dw.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=1, lty=2)
lines(exp(ww.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=1)
lines(exp(wd.pred.RatioL$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=1, lty=2)

plot(dryMass ~ numDaysField.unscaled, data = decomp, col = decomp.col, cex = 2, font.axis = 
       1, font.lab = 2, pch=pch, main="Title", ylab="Dry Mass",  xlab="Number of Days in Field", xlim=c(0,160), type="n")
lines(exp(dd.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="#0c7c59", lwd=2)
lines(exp(dd.pred.RatioH$mod.avg.pred) ~ newDays.unscaled, col="#0c7c59", lwd=2)



###run this 
pdf(file = "decomp.pdf", width = 12, height = 11.5, family = "serif")

plot(dryMass ~ numDaysField.unscaled, data = decomp, col = decomp.col, cex = 2,cex.lab=1.5, cex.axis=1.3, cex.main=1.8, cex.sub=1.5, font.axis = 
       1, font.lab = 2, pch=pch, main="Mass loss across 4 months", ylab="Dry Mass (g)",  xlab="Number of Days in Field", xlim=c(0,120), type="n",yaxt="n",family='serif')
p <- pretty(par("usr")[3:4])
l <- formatC(p,format="f", digits=2)
axis(2, at=p, labels=l, family="serif", cex.axis=1.2)


points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="wet" & Source=="wet"), col="DodgerBlue", pch=16, cex=1)
points(dryMass~numDaysField.unscaled, data=decomp, subset = c(SiteC=="dry" & Source=="wet"), col="pink", pch=17, cex=1)

points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="dry" & Source=="dry"), col="Brown", pch=2, cex=1)
points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="wet"& Source=="dry"), col="DarkGreen", pch=1, cex=1)

legend("bottomleft",legend=c("Dry litter in Dry plot", "Wet litter in Dry plot","Wet litter in Wet plot","Dry litter in Wet plot","Wet plot|Ratio:95th percentile","Wet plot|Ratio:5th percentile","Dry plot|Ratio:95th percentile","Dry plot|Ratio:5th percentile"),
       col=c("Brown", "pink","DodgerBlue","DarkGreen","blue","blue","red","red"), pch=c(2,17,16,1,NA,NA,NA,NA),lty=c(NA,NA,NA,NA,1,2,1,2),cex=1.2)


final.model=lme(LdryMass~
                  LinitialMass+numDaysField:SiteC+SiteC+Ratio + numDaysField, 
                random=list(~1|Site,~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	
                data=decomp, subset=complete.cases(decomp), na.action=na.exclude, method="ML")

new.wet.data.RatioH<-data.frame(SiteC="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.dry.data.RatioH<-data.frame(SiteC="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.wet.data.RatioL<-data.frame(SiteC="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.dry.data.RatioL<-data.frame(SiteC="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))

dry.RatioH.pred<-modavgPred(list(final.model), newdata=new.dry.data.RatioH)
wet.RatioH.pred<-modavgPred(list(final.model), newdata=new.wet.data.RatioH)
dry.RatioL.pred<-modavgPred(list(final.model), newdata=new.dry.data.RatioL)
wet.RatioL.pred<-modavgPred(list(final.model), newdata=new.wet.data.RatioL)


lines(exp(dry.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=1)

lines(exp(dry.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=2)

lines(exp(wet.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=1)

lines(exp(wet.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=2)
dev.off()

summary(final.model)
r.squaredGLMM(final.model)


#######
plot_model(final.model, type = "int", terms = c("numDaysField", "SiteC"),pred.type = c("fe", "re"))
?plot_model
#######

model1=lm(mpg~am*wt, data=dat)
summary(model1)

head(decomp)
decomp$predicted=predict(final.model)

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))

p=ggplot(dat, aes(x=wt, y=mpg, shape=am))+
  geom_point()+
  scale_shape_manual(values=c(1,16), name='Transmission', labels=c('Automatic','Manual'))+
  geom_line(aes(x = wt, y = predicted, linetype=am)) +
  scale_linetype_discrete(name='Transmission', labels=c('Automatic','Manual'))+
  labs(x = 'Vehicle Weight', y = 'Vehicle MPG')+
  apatheme
p


model<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:numDaysField
            , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)

model<-lme(fixed=log(dryMass) ~ numDaysField + Ratio + SiteC + SiteC:numDaysField
           , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)

model<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + SiteC + SiteC:numDaysField
           , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)

model<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + Ratio + SiteC + SiteC:numDaysField+ numDaysField:Ratio
           , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)

model<-lme(fixed=log(dryMass) ~ numDaysField + Ratio + SiteC + SiteC:numDaysField+ numDaysField:Ratio
           , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)

model<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + SiteC + Source + SiteC:numDaysField 
             , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)


model<-lme(fixed=log(dryMass) ~ numDaysField + SiteC + Source + SiteC:numDaysField 
           , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)

model<-lme(fixed=log(dryMass) ~ log(initialMass) + numDaysField + SiteC + Source + SiteC:numDaysField + numDaysField:Source
           , random=list(~ 1|Site, ~1|Batch), weights= varExp(form =~ as.vector(numDaysField) ),	data=decomp, na.action=na.exclude,method="ML")
r.squaredGLMM(model)


####extract values###
estDays <-(c(90,180,360)-mean(numDaysField.unscaled))/sd(numDaysField.unscaled)
newDays.unscaled<-(newDays*sd(numDaysField.unscaled)+mean(numDaysField.unscaled))


est.wet.data.RatioH<-data.frame(SiteC="wet",numDaysField=estDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
est.dry.data.RatioH<-data.frame(SiteC="dry",numDaysField=estDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
est.wet.data.RatioL<-data.frame(SiteC="wet",numDaysField=estDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
est.dry.data.RatioL<-data.frame(SiteC="dry",numDaysField=estDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))

est.dry.RatioH.pred<-modavgPred(list(final.model), newdata=est.dry.data.RatioH)
est.wet.RatioH.pred<-modavgPred(list(final.model), newdata=est.wet.data.RatioH)
est.dry.RatioL.pred<-modavgPred(list(final.model), newdata=est.dry.data.RatioL)
est.wet.RatioL.pred<-modavgPred(list(final.model), newdata=est.wet.data.RatioL)

exp(est.dry.RatioH.pred$mod.avg.pred)
exp(est.wet.RatioH.pred$mod.avg.pred)
exp(est.dry.RatioL.pred$mod.avg.pred)
exp(est.wet.RatioL.pred$mod.avg.pred)

summary(final.model)
###look here##


plot(dryMass ~ numDaysField.unscaled, data = decomp, col = decomp.col, cex = 2,cex.lab=1.5, cex.axis=1.3, cex.main=1.8, cex.sub=1.5, font.axis = 
       1, font.lab = 2, pch=pch, main="Mass loss across 4 months", ylab="Dry Mass (g)",  xlab="Number of Days in Field", xlim=c(0,360), type="n",yaxt="n",family='serif')
p <- pretty(par("usr")[3:4])
l <- formatC(p,format="f", digits=2)
axis(2, at=p, labels=l, family="serif", cex.axis=1.2)


points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="wet" & Source=="wet"), col="DodgerBlue", pch=16, cex=1)
points(dryMass~numDaysField.unscaled, data=decomp, subset = c(SiteC=="dry" & Source=="wet"), col="pink", pch=17, cex=1)

points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="dry" & Source=="dry"), col="Brown", pch=2, cex=1)
points(dryMass~numDaysField.unscaled, data=decomp, subset=c(SiteC=="wet"& Source=="dry"), col="DarkGreen", pch=1, cex=1)

legend("topright",legend=c("Dry litter in Dry plot", "Wet litter in Dry plot","Wet litter in Wet plot","Dry litter in Wet plot","Wet plot|Ratio:95th percentile","Wet plot|Ratio:5th percentile","Dry plot|Ratio:95th percentile","Dry plot|Ratio:5th percentile"),
       col=c("Brown", "pink","DodgerBlue","DarkGreen","blue","blue","red","red"), pch=c(2,17,16,1,NA,NA,NA,NA),lty=c(NA,NA,NA,NA,1,2,1,2),cex=1.2)

new.wet.data.RatioH<-data.frame(SiteC="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.dry.data.RatioH<-data.frame(SiteC="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.wet.data.RatioL<-data.frame(SiteC="wet",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.dry.data.RatioL<-data.frame(SiteC="dry",numDaysField=newDays,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))

dry.RatioH.pred<-modavgPred(list(final.model), newdata=new.dry.data.RatioH)
wet.RatioH.pred<-modavgPred(list(final.model), newdata=new.wet.data.RatioH)
dry.RatioL.pred<-modavgPred(list(final.model), newdata=new.dry.data.RatioL)
wet.RatioL.pred<-modavgPred(list(final.model), newdata=new.wet.data.RatioL)

lines(exp(dry.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=1)
lines(exp(dry.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=2)
lines(exp(wet.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=1)
lines(exp(wet.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=2)

#extrapolate
newDays360 <- seq(min(0-mean(numDaysField.unscaled))/sd(numDaysField.unscaled),max(360-mean(numDaysField.unscaled))/sd(numDaysField.unscaled), length.out = 1000)
newDays.unscaled<-(newDays360*sd(numDaysField.unscaled)+mean(numDaysField.unscaled))


new.wet.data.RatioH.360<-data.frame(SiteC="wet",numDaysField=newDays360,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.dry.data.RatioH.360<-data.frame(SiteC="dry",numDaysField=newDays360,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)+2*sd(decomp$Ratio))
new.wet.data.RatioL.360<-data.frame(SiteC="wet",numDaysField=newDays360,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))
new.dry.data.RatioL.360<-data.frame(SiteC="dry",numDaysField=newDays360,LinitialMass=mean(log(decomp$initialMass)), Ratio=mean(decomp$Ratio)-2*sd(decomp$Ratio))

dry.RatioH.pred<-modavgPred(list(final.model), newdata=new.dry.data.RatioH.360)
wet.RatioH.pred<-modavgPred(list(final.model), newdata=new.wet.data.RatioH.360)
dry.RatioL.pred<-modavgPred(list(final.model), newdata=new.dry.data.RatioL.360)
wet.RatioL.pred<-modavgPred(list(final.model), newdata=new.wet.data.RatioL.360)

lines(exp(dry.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=1)
lines(exp(dry.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="red", lwd=2, lty=2)
lines(exp(wet.RatioH.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=1)
lines(exp(wet.RatioL.pred$mod.avg.pred) ~ newDays.unscaled, col="blue", lwd=2, lty=2)


               