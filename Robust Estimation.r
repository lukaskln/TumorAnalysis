####### Robust Model #######

### Model with heteroscedastic Variance ###

fit.lme.ml.vexp.wogi <- lme(volume ~ daydiff + group:daydiff , random = ~ daydiff | mouse, data = dfData,
                                corr = corARMA(form= ~1|mouse, p = 2, q =0), # To model dependence. Unit root in MA for ARMA(2,1)
                                weights = varExp(form=~fitted(.)|v3SplitGroups), # Three groups can be changed to two as well
                                method = "ML", na.action =na.pass)


summary(fit.lme.ml.vexp.wogi)

pred.ml.vexp.wogi = predict(fit.lme.ml.vexp.wogi)

res.ml.vexp.wogi.pea = resid(fit.lme.ml.vexp.wogi, type = "p")

res.ml.vexp.wogi = resid(fit.lme.ml.vexp.wogi)

p13 =
ggplot(data = dfData,aes(y = res.ml.vexp.wogi.pea, x = pred.ml.vexp.wogi)) +
  geom_point(aes(color = group)) +
  scale_color_manual(values=vColors[1:2]) +
  geom_line(y = 0) +
  theme_minimal() + 
  xlab("Fitted by ML and with VarExp & CorARMA modeled") + 
  ylab("Pearson Residuals")

par(mfrow=c(1,1))

qqnorm(resid(fit.lme.ml.vexp.wogi, type ="p"),main="Pearson Residuals",pch=16,bty='n', col = vColors[3])
qqline(resid(fit.lme.ml.vexp.wogi, type = "p"),lwd=2,col= vColors[2]) 

ggplot(data = dfData,aes(y = res.ml.vexp.wogi.pea, x = v3SplitGroups)) +
  geom_boxplot(fill = vColors[1:3], alpha = 0.6, color= vColors[1:3]) +
  theme_minimal() + 
  xlab("Three Groups") + 
  ylab("Pearson Residuals fitted by ML and with VarExp & CorARMA modeled")

anova(fit.lme.ml.vexp.wogi, fit.lme.ml.wogi)


### Model with forced uncorrelation between random intercept and slope ###

summary(
  lme(volume ~ daydiff + group:daydiff, 
  random = list(mouse=pdIdent(~daydiff)), 
  data = dfData,
  method = "ML")
  ) # Not good


### Robust model with HUber Weights ###

fit.rob.ml.wogi <- rlmer(volume ~ daydiff + group:daydiff + (daydiff | mouse), data = dfData, 
                         rho.sigma.b=psi2propII(smoothPsi, k=1.345), # Higher k more efficiency, lower more robust
                         rho.sigma.e=psi2propII(smoothPsi, k=1.345)
                         )

summary(fit.rob.ml.wogi)

Weights = summary(fit.rob.ml.wogi)[9]$wgt.e

pred.rob.ml.wogi = predict(fit.rob.ml.wogi) 

res.rob.ml.wogi = residuals(fit.rob.ml.wogi, type = "response", scaled = TRUE)

p14 = 
ggplot(data = dfData,aes(y = res.rob.ml.wogi, x = pred.rob.ml.wogi, color = Weights)) +
  geom_point() +
  scale_color_gradient(high = vColors[2], low = vColors[1]) +
  geom_line(y = 0) +
  theme_minimal() + 
  xlab("Fitted by ML and Huber Loss function") + 
  ylab("Scaled Residuals w/o weights applied")

qqnorm(resid(fit.rob.ml.wogi),main="Residuals",pch=16,bty='n', col = vColors[3])
qqline(resid(fit.rob.ml.wogi),lwd=2,col= vColors[2]) 

compare(fit.lme4.ml.wogi, fit.rob.ml.wogi)

grid.arrange(p13,p14, nrow = 1)

plot.rlmerMod(fit.rob.ml.wogi)
