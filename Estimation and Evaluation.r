#######  Model Estimation #######

### Write model equation for full and maybe some restricted models to compare. Watch out that model matches research hypothesis! ###

# Full
fit.lme.ml.full <- lme(volume ~ group*daydiff, random = ~ daydiff | mouse, data = dfData, method = "ML")
summary(fit.lme.ml.full)

fit.lme.reml.full <- lme(volume ~ group*daydiff, random = ~ daydiff | mouse, data = dfData, method = "REML") # REML Estimation
summary(fit.lme.reml.full)

# Without Random Slope
fit.lme.ml.wors <- lme(volume ~ group*daydiff, random = ~ 1 | mouse, data = dfData, method = "ML")
summary(fit.lme.ml.wors)

# Without Group intercept
fit.lme.ml.wogi <- lme(volume ~ daydiff + group:daydiff , random = ~ daydiff | mouse, data = dfData, method = "ML")
summary(fit.lme.ml.wogi)

fit.lme.reml.wogi <- lme(volume ~ daydiff + group:daydiff , random = ~ daydiff | mouse, data = dfData, method = "REML")
summary(fit.lme.reml.wogi)

# Without Group intercept and without Random Slope
fit.lme.ml.wogirs <- lme(volume ~ daydiff + group:daydiff, random = ~ 1 | mouse, data = dfData, method = "ML")
summary(fit.lme.ml.wogirs)

# Independent Random intercept & slope
fit.lme.ml.wogi.random.independent = lme(volume ~ daydiff + group:daydiff, data = dfData, random = list(mouse = pdDiag(~daydiff)), method = "ML")
summary(fit.lme.ml.wogi.random.independent)


####### Model Validation #######

### Check model performance based on IC and CIs ###

anova(fit.lme.ml.full, fit.lme.ml.wogi)

head(model.matrix(lmer(volume ~ daydiff*group + (daydiff | mouse), data = dfData, REML = FALSE)),10) # Plug in different Models

intervals(fit.lme.ml.full)
intervals(fit.lme.ml.wogi)

anova(fit.lme.ml.full, type="sequential" )

round(r.squaredGLMM(fit.lme.ml.wogi),4)
round(r.squaredGLMM(fit.lme.ml.full),4)


### Predictions ###

pred.ml.wogi = predict(fit.lme.ml.wogi,level=1)

res.ml.wogi = resid(fit.lme.ml.wogi)

pred.reml.wogi = predict(fit.lme.reml.wogi, level = 1)

res.reml.wogi = resid(fit.lme.reml.wogi)

### Check prediction (errors) for: homoscedasticity, outliers, normality and linearity (Correlation between random effects?). ###
### Also within factors (Scatterplot, qqplot). ###

p8=
ggplot(data = dfData,aes(y = res.ml.wogi, x = pred.ml.wogi)) +
  geom_point(color = vColors[1]) +
  geom_line(y = 0) +
  theme_minimal() + 
  xlab("Fitted by ML") + 
  ylab("Residuals")

p9=
  ggplot(data = dfData,aes(y = res.reml.wogi, x = pred.ml.wogi)) +
  geom_point(color = vColors[2]) +
  geom_line(y = 0) +
  theme_minimal() + 
  xlab("Fitted by REML") + 
  ylab("Residuals")

p8.1=
  ggplot(data = dfData,aes(y = res.ml.wogi, x = pred.ml.wogi)) +
  geom_point(aes(color = group)) +
  scale_color_manual(values=vColors[1:2]) +
  geom_line(y = 0) +
  theme_minimal() + 
  xlab("Fitted by ML") + 
  ylab("Residuals")

p9.1=
  ggplot(data = dfData,aes(y = res.reml.wogi, x = pred.ml.wogi)) +
  geom_point(aes(color = v3SplitGroups)) +
  scale_color_manual(values=vColors[1:3]) +
  geom_line(y = 0) +
  theme_minimal() + 
  xlab("Fitted by ML") + 
  ylab("Residuals")

grid.arrange(p8,p9,p8.1, p9.1, nrow = 2)

p10 =
  ggplot(data = dfData,aes(y = res.ml.wogi, x = group)) +
  geom_boxplot(fill = vColors[1:2], alpha = 0.6, color= vColors[1:2]) +
  theme_minimal() + 
  xlab("Two Groups") + 
  ylab("Residuals fitted by ML")

p10.1 =
  ggplot(data = dfData,aes(y = res.ml.wogi, x = v3SplitGroups)) +
  geom_boxplot(fill = vColors[1:3], alpha = 0.6, color= vColors[1:3]) +
  theme_minimal() + 
  xlab("Three Groups") + 
  ylab("Residuals fitted by ML")

p11 =
  ggplot(data = dfData,aes(y = res.ml.wogi, x = as.factor(daydiff))) +
  geom_boxplot(fill = vColors[3], alpha = 0.6, color= vColors[2]) +
  theme_minimal() + 
  xlab("Days Passed") + 
  ylab("Residuals fitted by ML")

p12 =
  ggplot(data = dfData,aes(y = res.ml.wogi, x = mouse)) +
  geom_boxplot(fill = vColors[3], alpha = 0.6, color= vColors[1]) +
  theme_minimal() + 
  xlab("Mouse") + 
  ylab("Residuals fitted by ML")

grid.arrange(p10,p10.1, p11,p12, layout_matrix = rbind(c(1,2),c(3,4)))

### Compute Block-Correlation Matrix ###

Omega.matrix = function(model) {
  if(class(model)[1]=="lmerModLmerTest"|class(model)[1]=="lmerMod"){
    sigma2.epsilon = sigma(model)^2
    Psi.star       = crossprod(getME(model,"Lambdat"))*sigma2.epsilon
    Z              = getME(model,"Z")
    Omega          = Z %*% Psi.star %*% t(Z) + sigma2.epsilon* Diagonal(nrow(Z))
    Omega
  }else{
    warning("Function only works on outcput of function lmer()\n")
  }
}

Omega.plot = function(Omega,legend=TRUE,axes=TRUE){
  corw = cov2cor(Omega)
  if(any(corw<0)){
    colw=c(gray(1),rainbow(197)[197:100],gray(.9),rainbow(197)[99:1],gray(0))
    zlim=c(-1,1)
  }else{
    colw=c(gray(.9),rainbow(98)[98:1],gray(0))
    zlim=c(0,1)
  }
  image(z=as.matrix(corw[nrow(corw):1,]),zlim=zlim,axes=FALSE,col=colw)
  if(legend){    
    valw = as.numeric(names(table(as.matrix(corw))))  
    posw = round(valw*length(colw))
    posw[posw==0] = 1
    posw[posw>length(colw)] = length(colw)    
    legend("topright",ncol=1,legend=format(round(valw,4)),
           col=colw[posw],pch=15,bg="light gray",
           title="Values",box.lwd=NA)
  }
  if(axes){
    axis(3,at=seq(0,1,length=nrow(Omega)),labels=FALSE)
    axis(2,at=seq(0,1,length=nrow(Omega)),labels=FALSE)
    axis(2,at=c(1,0),c(1,nrow(Omega)),las=2)
    axis(3,at=c(0,1),c(1,nrow(Omega)),las=1)
  }
}

fit.lme4.ml.wogi <- lmer(volume ~ daydiff + group:daydiff + (daydiff| mouse), data = dfData, REML = FALSE) #lme4 object needed


sigma2.epsilon = sigma(fit.lme4.ml.wogi)^2
Psi.star       = crossprod(getME(fit.lme4.ml.wogi,"Lambdat"))*sigma2.epsilon
Z              = getME(fit.lme4.ml.wogi,"Z")
Omega          = Z %*% Psi.star %*% t(Z) + sigma2.epsilon* Diagonal(nrow(Z))

D = sqrt(diag(diag(Omega)))
DInv = solve(D)
CorOmega = DInv %*% Omega %*% DInv


Omega = Omega.matrix(fit.lme4.ml.wogi)

par(mfrow= c(1,2))
Omega.plot(Omega)
Omega.plot(Omega[-c(1:214),-c(1:214)], legend = FALSE)

### (Partial) Autocorrelation Function ###
par(mfrow= c(1,2))
acf(residuals(fit.lme.ml.wogi))
pacf(residuals(fit.lme.ml.wogi))

### Check normal distribution assumption for residuals and random effects (qqplot) ###

par(mfrow=c(1,3))

qqnorm(ranef(fit.lme.ml.wogi)[,1],main="Random intercept",pch=16,bty='n', col = vColors[3])
qqline(ranef(fit.lme.ml.wogi)[,1],lwd=2,col= vColors[1])

qqnorm(ranef(fit.lme.ml.wogi)[,2],main="Random slope",pch=16,bty='n', col = vColors[3])
qqline(ranef(fit.lme.ml.wogi)[,2],lwd=2,col= vColors[1])

qqnorm(resid(fit.lme.ml.wogi),main="Residuals",pch=16,bty='n', col = vColors[3])
qqline(resid(fit.lme.ml.wogi),lwd=2,col= vColors[2])

ggplot(data = ranef(fit.lme.ml.wogi),aes(y = ranef(fit.lme.ml.wogi)[,2], x = ranef(fit.lme.ml.wogi)[,1])) +
  geom_point(color= vColors[3]) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, color = vColors[2]) +
  theme_minimal() + 
  xlab("Random Intercept") + 
  ylab("Random Slope")
