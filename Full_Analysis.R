library(lme4)
library(lmerTest)
library(nlme)
library(scales)
library(readr)
library(forecast)
library(robustlmm)
library(MuMIn)
library(mvtnorm)

# 1. Data Transformation

TumourSize<- read.table("TumourSize.txt",header=T,sep=",") # Read Dataset 

TumourSize = TumourSize[-1] # Remove X column

Vol = TumourSize$volume # Non transformed Volume
ThreeGroups = TumourSize$group # Not combined groups

TumourSize$volume = log(TumourSize$volume)  # apply log to volume
TumourSize$date = as.Date(TumourSize$date) # Change from factor to date
TumourSize["dayspassed"] = NA # New column with days difference from starting date

for(i in levels(TumourSize$mouse)){
  
  TumourSize$dayspassed [TumourSize$mouse == i] = 
    difftime( TumourSize$date[TumourSize$mouse == i], 
              min(TumourSize$date[TumourSize$mouse == i]), unit = "days")
  
}

TumourSize$group = as.character(TumourSize$group)
TumourSize$group[which(TumourSize$group == "Control I" | TumourSize$group == "Control II")] = "Control" # Set both groups equal
TumourSize$group = as.factor(TumourSize$group)
contrasts(TumourSize$group) = contr.treatment(2, base = 1) # "Control" is now baseline group

summary(TumourSize) # Unbalanced distribution of observations inside mouse and group 


# Descriptive Analysis

# QQ Plots before and after log transformation of Volume
par(mfrow=c(1,2))

qqnorm(Vol,pch=16,bty='n',main='Volume', col = "gray50", ylab = "Quantiles of mm³")
qqline(Vol,lwd=2,col= "blue")

qqnorm(TumourSize$volume,pch=16,bty='n',main='Log-Volume', col = "gray50", ylab = "Quantiles of mm³")
qqline(TumourSize$volume,lwd=2,col= "green")

# Proof of linear Relationship
#Boxplot
par(mfrow=c(1,4))

boxplot(volume~ThreeGroups,main="group",data=TumourSize,cex.axis=1,col="light gray")
boxplot(volume~group,main="group",data=TumourSize,cex.axis=1,col="light gray")
boxplot(volume~mouse,main="mouse",data=TumourSize,cex.axis=1,col="light gray")
boxplot(volume~dayspassed,main="date",data=TumourSize,cex.axis=1,col="light gray")

summary(aov(volume ~ group, data = TumourSize))

summary(aov(volume ~ ThreeGroups, data = TumourSize))

summary(aov(volume ~ ThreeGroups[ThreeGroups != "Treatment"], data = TumourSize[TumourSize$group != "Treatment",])) # Anova between the two controll groups


par(mfrow=c(1,1))
plot(x = dayspassed, y = volume, col = ifelse(group == "Treatment", "red", "blue"), pch = 16)
abline(lm(volume ~ dayspassed, data = TumourSize[TumourSize$group != "Treatment",]), col = "blue", lwd = 2) # Control group Line
abline(lm(volume ~ dayspassed, data = TumourSize[TumourSize$group == "Treatment",]), col = "red", lwd = 2) # Treatment Group Line

plot(groupedData(volume~ dayspassed|mouse,data= TumourSize,order.groups=F),inner = TumourSize$group, ylab = "log of Volume")
plot(groupedData(volume~ dayspassed|group,data= TumourSize,order.groups=F), ylab = "log of Volume")
plot(groupedData(volume~ dayspassed|group,data= TumourSize,order.groups=F),inner = TumourSize$mouse, ylab = "log of Volume")

#interaction of tow independed objects with the depended one
interaction.plot(group, dayspassed,volume, legend = FALSE)



#GroupedData: groupedData(Y~X1|X2)lustration of Y vs X1(fixed) for each X2(random)
TS.plot1<-groupedData(volume~group|dayspassed,data=TumourSize,order.groups=F)
plot(TS.plot1)
plot(TS.plot1,outer=~group)
TS.plot2<-groupedData(volume~date|mouse,data=TumourSize,order.groups=F)
plot(TS.plot2)
plot(TS.plot2,outer=~group)
TS.plot4<-groupedData(volume~mouse|date,data=TumourSize,order.groups=F)
plot(TS.plot4)
plot(TS.plot4,outer=~group)

#NO Plot with Jitter (Add Noise to Numbers) SINCE NO numeric indipedent variable



#Column rank design matrix for Group (Check and change levels and contrasts of variable Group)
#TumourSize$group = factor(TumourSize$group, levels=c("Treatment","Control I","Control II"))
contrasts(TumourSize$group)
model.matrix(~TumourSize$group,group)
#contrasts(TumourSize$group) <- contr.treatment(3, base = 1)

# 3. Models

TumourSize.lme1.ml <- lme(fixed = volume~ group + dayspassed + group:dayspassed,
                       data=TumourSize,random= ~dayspassed|mouse,
                       method = "ML")
summary(TumourSize.lme1.ml) # interaction with the fixed effects

TumourSize.lme1.reml <- lme(fixed = volume~ group + dayspassed + group:dayspassed,
                          data=TumourSize,random= ~dayspassed|mouse,
                          method = "REML")
summary(TumourSize.lme1.reml)

# Y = beta0 + beta1*G + beta2*DaysPassed +beta3 * G * DaysPassed + RM1 + RM2 * DaysPassed + Error

TumourSize.lme2 <- lme(fixed = volume~  dayspassed + group:dayspassed,data=TumourSize,random= ~dayspassed|mouse, method = "ML")
summary(TumourSize.lme2)

anova(TumourSize.lme2, type="sequential" )

# Y = beta0 + beta2*DaysPassed +beta3 * G * DaysPassed + RM1 + RM2 * DaysPassed + Error

TumourSize.lme3 <- lme(fixed = volume~  dayspassed + group:dayspassed,data=TumourSize,
                       random= list(mouse = pdDiag(~dayspassed)),
                       method = "ML")
summary(TumourSize.lme3)

# Whats the difference with this: ?? TumourSize.lme3 <- lme(fixed=volume ~ group*date, random = ~ group | mouse, data = TumourSize)



# Model Validation

Predictions = predict(TumourSize.lme2)

Residuals = residuals(TumourSize.lme2, type = "p")

plot(TumourSize.lme2)
plot(TumourSize.lme2, resid(., type = "p") ~ fitted(.) | group)
plot(TumourSize.lme2, resid(., type = "p") ~ fitted(.) | mouse)

plot(ACF(TumourSize.lme2, plot = TRUE), alpha = 0.05)

boxplot(Residuals ~ groups)
boxplot(Residuals ~ ThreeGroups)
boxplot(Residuals ~ mouse)

qqnorm(Residuals,pch=16,bty='n',main='Residuals', col = "gray50", ylab = "Quantiles of mm³")
qqline(Residuals,lwd=2,col= "blue")


# PREDICTIONS
# Estimation of the random effects (not std. dev!)
estim.TS.ranef<-random.effects(TumourSize.lme2)
head(estim.TS.ranef)
## Predictions
pred01 <- predict(TumourSize.lme2,level=0:1)
head(pred01)

# Checking normality assumptions of the random effects
qqnorm(TumourSize.lme1,~ranef(.),id=0.05)
plot(TumourSize.lme1)

qqnorm(TumourSize.lme2,~ranef(.),id=0.05)
qqnorm(TumourSize.lme3,~ranef(.),id=0.05)

qqnorm(estim.TS.ranef[,1],datax=T,pch=16,bty='n',main='(intercept)')
qqline(estim.TS.ranef[,1],datax=T,lwd=1.5,col='red')


# MODEL CHECKING
# Fit random intercept and random slope model with REML
fit.lmer.reml <- lmer(volume ~ group*date + (1 | mouse), data = TumourSize)
summary(fit.lmer.reml)

#Models1
Model1 = lm(volume~mouse+date+group,data = TumourSize)
plot(Model1)
summary(Model1)

ModelLog1 = lm(log(volume)~mouse+date+group,data = TumourSize)
plot(ModelLog1)
summary(ModelLog1)

Model2 = lm(volume~mouse+date+group+mouse*group,data = TumourSize)
plot(Model2)
summary(Model2)

ModelLog2 = lm(log(volume)~mouse+date+group+mouse*group,data = TumourSize)
plot(ModelLog2)
summary(ModelLog2)

Model3 <- lme(volume~group+date,data=TumourSize,random=~1|mouse,method="ML")
plot(Model3)
summary(Model3)
qqnorm(Model3,~ranef(.),id=0.05)
pairs(Model3,~ranef(.))


# Robust model estimation 
TumourSize.lme4 =
lme(volume ~ dayspassed + group:dayspassed, 
    random = ~ dayspassed | mouse, data = TumourSize,
    weights = varPower(form=~fitted(.)|ThreeGroups), 
    method = "ML", na.action =na.pass)

summary(TumourSize.lme4)

Residuals4 = residuals(TumourSize.lme4, type = "p")

plot(TumourSize.lme4)
boxplot(Residuals4 ~ ThreeGroups)
