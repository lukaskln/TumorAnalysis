####### Setup #######

library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(lme4)
library(lmerTest)
library(nlme)
library(scales)
library(readr)
library(forecast)
library(robustlmm)
library(MuMIn)
library(mvtnorm)

vColors = c("#eb0087","#554696","#d2d2d2","#08919E")


####### Data Transformation #######
  
dfData = read.delim("TumourSize.txt", sep = ",")

dfData = dfData[-1] # Remove X column
  
dfData$volume = dfData$volume %>% log() # apply log to volume

dfData$date = dfData$date %>% as.Date() # Change from factor to date

dfData["daydiff"] = NA # New column with days difference from starting date

for(i in levels(dfData$mouse)){
  
  dfData$daydiff[dfData$mouse == i] = difftime( dfData$date[dfData$mouse == i], min(dfData$date[dfData$mouse == i]), unit = "days")
  
}

### Section for Plots before combining groups ###

p0 = 
ggplot(data = dfData,aes(y = volume, x = group)) +
  geom_boxplot(fill = vColors[1:3], alpha = 0.6, color= vColors[1:3]) +
  theme_minimal() + 
  xlab("Group") + 
  ylab("mm³ Volume (log)") # Plot for Volume distr. inside of groups

lAnova = summary(aov(volume ~ group, data = dfData))

lAnovaControl = summary(aov(volume ~ group, data = dfData[dfData$group != "Treatment",]))

###

dfData$group = as.character(dfData$group)

v3SplitGroups = dfData$group

dfData$group[which(dfData$group == "Control I" | dfData$group == "Control II")] = "Control" # Set both groups equal

dfData$group = as.factor(dfData$group)

contrasts(dfData$group) = contr.treatment(2, base = 1) # "Control" is now baseline group

summary(dfData)

####### Descriptive Analysis #######

### Linearity of the relation between the covariates and the response (Correlation) ###

p1 =
  ggplot(data = dfData,aes(y = volume, x = daydiff, color = v3SplitGroups)) +
  geom_point() +
  scale_color_manual(values=vColors[1:4])+
  geom_smooth(aes(color = v3SplitGroups), method='lm',formula=y~x, se=FALSE) +
  theme_minimal() + 
  xlab("Passed Days") + 
  ylab("mm³ Volume (log)")

p2 =
  ggplot(data = dfData,aes(y = volume, x = daydiff, color = group)) +
  geom_point() +
  scale_color_manual(values=vColors[1:2])+
  geom_smooth(aes(color = group), method='lm',formula=y~x, se=FALSE) +
  theme_minimal() + 
  xlab("Passed Days") + 
  ylab("mm³ Volume (log)")


grid.arrange(p1,p2, nrow = 1)

cor(dfData[c("volume", "daydiff")])

### Equal variance of the responses around the different factor's levels (Boxplot) (Anova) ###

lAnova # Anova for 3 groups
lAnovaControl # Anova only between Control Group 1 and 2
summary(aov(volume ~ group, data = dfData)) #Anova for 2 Groups
summary(aov(volume ~ mouse, data = dfData)) #Anova for Mice

p3 =
  ggplot(data = dfData,aes(y = volume, x = group)) +
  geom_boxplot(fill = vColors[1:2], alpha = 0.6, color= vColors[1:2]) +
  theme_minimal() + 
  xlab("Group") + 
  ylab("mm³ Volume (log)")

p4=
  ggplot(data = dfData,aes(y = volume, x = mouse)) +
  geom_boxplot(fill = vColors[3], alpha = 0.6, color= vColors[2]) +
  theme_minimal() + 
  xlab("Mouse") + 
  ylab("mm³ Volume (log)")

grid.arrange(p0, p3,p4, layout_matrix = rbind(c(1,2),c(3,3)))

### Presence or absence of interactions between the different "explanatory variables" (Interaction plot) ###

par(mfrow=c(1,2))
interaction.plot(dfData$group, dfData$daydiff, response = dfData$volume, legend = FALSE)
interaction.plot(dfData$mouse, dfData$daydiff, response = dfData$volume, legend = FALSE)

### Correlation Analysis for Random effects ###

plot(groupedData(volume~ daydiff|group,data=dfData,order.groups=F),inner = dfData$mouse, ylab = "log of Volume", key = FALSE)

### Outliers (Boxplot or Scatterplot) ###

p5 =
  ggplot(data = dfData,aes(y = volume, x = group)) +
  geom_jitter(fill = vColors[1], alpha = 0.6, color= vColors[1]) +
  theme_minimal() + 
  xlab("Group") + 
  ylab("mm³ Volume (log)")

p6 =
  ggplot(data = dfData,aes(y = volume, x = as.factor(daydiff))) +
  geom_boxplot(fill = vColors[2], alpha = 0.6, color= vColors[2]) +
  theme_minimal() + 
  xlab("Passed Days") + 
  ylab("mm³ Volume (log)")

p7 =
  ggplot(data = dfData,aes(y = volume, x = mouse)) +
  geom_boxplot(fill = vColors[3], alpha = 0.6, color= vColors[1]) +
  theme_minimal() + 
  xlab("Mouse") + 
  ylab("mm³ Volume (log)")

grid.arrange(p5, p6,p7, layout_matrix = rbind(c(1,2),c(3,3)))

### Normality (Given a Group or Marginal) (qqplot within groups) ###

par(mfrow=c(1,4))

qqnorm(exp(dfData$volume),pch=16,bty='n',main='Volume', col = vColors[3], ylab = "Quantiles of mm³")
qqline(exp(dfData$volume),lwd=2,col= vColors[1])

qqnorm(dfData$volume,pch=16,bty='n',main='Log Volume', col = vColors[3], ylab = "Quantiles of log mm³")
qqline(dfData$volume,lwd=2,col= vColors[1])

qqnorm(dfData$volume[dfData$group == "Control"],pch=16,bty='n',main='Log Volume within Control', col = vColors[3], ylab = "Quantiles of log mm³ for Control Group")
qqline(dfData$volume[dfData$group == "Control"],lwd=2,col= vColors[2])

qqnorm(dfData$volume[dfData$group == "Treatment"],pch=16,bty='n',main='Log Volume within Treatment', col = vColors[3], ylab = "Quantiles of log mm³ for Treatment Group")
qqline(dfData$volume[dfData$group == "Treatment"],lwd=2,col= vColors[2])
