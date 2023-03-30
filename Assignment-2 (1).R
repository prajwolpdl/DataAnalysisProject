install.packages('NHANES') ; install.packages('tidyverse')
library(tidyverse)
library(NHANES)
library(car)
library(MASS) 
library(glmnet)
library(rms)
small.nhanes <- na.omit(NHANES[NHANES$SurveyYr=="2011_12"
                               & NHANES$Age > 17,c(1,3,4,8:11,13,17,20,21,25,46,50,51,52,61)])
small.nhanes <- as.data.frame(small.nhanes %>%
                                group_by(ID) %>% filter(row_number()==1) )
nrow(small.nhanes)

set.seed(1002656486)
train <- small.nhanes[sample(seq_len(nrow(small.nhanes)), size = 400),]
nrow(train)
length(which(small.nhanes$ID %in% train$ID))
test <- small.nhanes[!small.nhanes$ID %in% train$ID,]
nrow(test)

count_smokenow <- table(small.nhanes$SmokeNow)
count = as.data.frame(count_smokenow)
count

full_model <- lm(BPSysAve ~ Gender + Age + Race3 + Education + MaritalStatus + HHIncome + Poverty + Weight + Height + BMI + Depressed + SleepTrouble + PhysActive + SmokeNow, data = train )
summary(full_model)
anova(full_model)
vif(full_model)

## BPSys vs SmokeNow ##

boxplot(train$BPSysAve~train$SmokeNow, main = "Boxplot for blood pressure vs smoking status", xlab = "Smoking Status", ylab ="Systolic BP levels")
hist(train$BPSysAve, breaks= "Scott", prob=TRUE, xlab="Systolic Blood Pressure levels", ylab="Density", main= "Histogram for combined systolic Blood Pressure reading")
x1 <- data.frame(train$BPSysAve, train$SmokeNow)
x1

t.test(train$BPSysAve~train$SmokeNow, data = x1, var.equal = TRUE)


## Step-wise Regression ##

##AIC##
sel.var.aic.mod <- step(full_model, trace = 0, k = 2, direction = "both") 
sel.var.aic<-attr(terms(sel.var.aic.mod), "term.labels")   
sel.var.aic

##BIC##

n <- nrow(train)
sel.var.bic.mod <- step(full_model, trace = 0, k = log(n), direction = "both") 
sel.var.bic<-attr(terms(sel.var.bic.mod), "term.labels")   
sel.var.bic

## diagnostics checking after stepwise regression ##

vif(sel.var.aic.mod)
vif(sel.var.bic.mod)

## The hat values ###

h.aic <- hatvalues(sel.var.aic.mod)
thresh <- 2 * (dim(model.matrix(sel.var.aic.mod))[2])/nrow(train)
w.aic <- which(h.aic > thresh)
w.aic

h.bic <- hatvalues(sel.var.bic.mod)
thresh <- 2 * (dim(model.matrix(sel.var.bic.mod))[2])/nrow(train)
w.bic <- which(h.bic > thresh)
w.bic

## Standardized residuals ##

aic_mod.stres <- rstudent(sel.var.aic.mod)

fitted <- predict(sel.var.aic.mod)
qqnorm(aic_mod.stres)
plot(aic_mod.stres ~ fitted, type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

bic_mod.stres <- rstudent(sel.var.bic.mod)

fitted <- predict(sel.var.bic.mod)
qqnorm(bic_mod.stres)
plot(bic_mod.stres ~ fitted, type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

## Cook's distance##
D.aic <- cooks.distance(sel.var.aic.mod)
D.aic
which(D.aic > qf(0.5, 5, 168-5))

D.bic <- cooks.distance(sel.var.bic.mod)
D.bic
which(D.bic > qf(0.5, 5, 168-5))

## DFFITS ##
dfits.aic <- dffits(sel.var.aic.mod)
which(abs(dfits.aic) > 2*sqrt(5/168))

dfits.bic <- dffits(sel.var.bic.mod)
which(abs(dfits.bic) > 2*sqrt(5/168))

## DFBETAS ##
dfb.aic <- dfbetas(sel.var.aic.mod)
which(abs(dfb.aic[,2]) > 2/sqrt(168))

dfb.bic <- dfbetas(sel.var.bic.mod)
which(abs(dfb.bic[,2]) > 2/sqrt(168))



### Cross Validation and prediction performance of AIC based selection ###
ols.aic <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.aic, "BPSysAve"))], 
               x=T, y=T, model = T)

## 10 fold cross validation ##    
aic.cross <- calibrate(ols.aic, method = "crossvalidation", B = 10)
## Calibration plot ##
pdf("aic_cross.pdf", height = 8, width = 16)
plot(aic.cross, las = 1, xlab = "Predicted Blood Pressure reading", main = "Cross-Validation calibration with AIC")
dev.off()
## Test Error ##
pred.aic <- predict(ols.aic, newdata = test[,which(colnames(train) %in% c(sel.var.aic, "BPSysAve"))])
## Prediction error ##
pred.error.AIC <- mean((test$BPSysAve - pred.aic)^2)


### Cross Validation and prediction performance of BIC based selection ###
ols.bic <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.bic, "BPSysAve"))], 
               x=T, y=T, model = T)


## 10 fold cross validation ##    
bic.cross <- calibrate(ols.bic, method = "crossvalidation", B = 10)
## Calibration plot ##
pdf("bic_cross.pdf", height = 8, width = 16)
plot(bic.cross, las = 1, xlab = "Predicted Blood Pressure reading", main = "Cross-Validation calibration with BIC")
dev.off()
sel.var.bic<-attr(terms(ols.bic), "term.labels") 
sel.var.bic


## Test Error ##
pred.bic <- predict(ols.bic, newdata = test[,which(colnames(train) %in% c(sel.var.bic, "BPSysAve"))])
## Prediction error ##
pred.error.BIC <- mean((test$BPSysAve - pred.bic)^2)

##LASSO cross validation ##
set.seed(1004993297)
factors <- model.matrix(BPSysAve ~ Gender + Race3 + Education + MaritalStatus + HHIncome + Depressed + SleepTrouble + PhysActive + SmokeNow + Age, data = train)[, -1]
X <- as.matrix(data.frame(train[,8:11], train$SleepHrsNight, factors))

model.lassocv <- cv.glmnet(x = X, y = train$BPSysAve, standardize = T, alpha = 1)
best.lassocv <- model.lassocv$lambda.min
co<-coef(model.lassocv, s = "lambda.1se")

## threshold for variable selection ##

thresh <- 0.00
# select variables #
inds<-which(abs(co) > thresh )
variables<-row.names(co)[inds]
sel.var.lasso<-variables[!(variables %in% '(Intercept)')]
sel.var.lasso

ols.lasso <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.lasso, "BPSysAve"))], 
                 x=T, y=T, model = T)

## 10 fold cross validation ##    
lasso.cross <- calibrate(ols.lasso, method = "crossvalidation", B = 10)
## Calibration plot ##
pdf("lasso_cross.pdf", height = 8, width = 16)
plot(lasso.cross, las = 1, xlab = "Predicted LPSA", main = "Cross-Validation calibration with LASSO")
dev.off()

## Test Error ##
pred.lasso <- predict(ols.lasso, newdata = test[,which(colnames(train) %in% c(sel.var.lasso, "BPSysAve"))])
## Prediction error ##
pred.error.lasso <- mean((test$BPSysAve - pred.lasso)^2)

transformed <- lm(cbind(train$BPSysAve ,X) ~ 1)
bc <- powerTransform(transformed)
summary(bc)

## The hat values ###

h <- hatvalues(full_model)
thresh <- 2 * (dim(model.matrix(full_model))[2])/nrow(train)
w <- which(h > thresh)
w

## Standardized residuals ##


full_model.stres <- rstudent(full_model)

fitted <- predict(full_model)
qqnorm(full_model.stres, main = "Q-Q plot of final model", ylab = "Standardized residuals")
plot(full_model.stres ~ fitted, type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

## Cook's distance##
D <- cooks.distance(full_model)
D
which(D > qf(0.5, 5, 168-5))

## DFFITS ##
dfits <- dffits(full_model)
which(abs(dfits) > 2*sqrt(5/168))

## DFBETAS ##
dfb <- dfbetas(full_model)
which(abs(dfb[,2]) > 2/sqrt(168))


## Final model ##

summary(sel.var.bic.mod)

final.mod <- lm(BPSysAve~ Age+Poverty+Weight+SleepTrouble+SmokeNow, data= train)
summary(final.mod)
vif(final.mod)

par(mfrow = c(1,2))

final_mod.stres <- rstudent(final.mod)

fitted.final <- predict(final.mod)
qqnorm(final_mod.stres)
plot(final_mod.stres ~ fitted.final, main= "Residual plot for the final model", type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

mult <- boxcox(train$BPSysAve ~ train$Age + train$Poverty + train$Weight + train$SleepTrouble + train$SmokeNow)
lambda <- mult$x[which.max(mult$y)]
transform_mod <- lm(((BPSysAve^lambda-1)/lambda)~ Age + Poverty + Weight + SleepTrouble + SmokeNow, data = train)

transform_mod.stres <- rstudent(transform_mod)

fitted.transform <- predict(transform_mod)
qqnorm(transform_mod.stres)
plot(transform_mod.stres ~ fitted.transform, main="Residual plot for the transformed model", type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

sel.var.final<-attr(terms(final.mod), "term.labels")
sel.var.final

sel.var.transform <- attr(terms(transform_mod), "term labels")
sel.var.transform

summary(transform_mod)
vif(transform_mod)

## cross-validation ##

newX <- model.matrix(BPSysAve ~ Weight+ Poverty + SmokeNow + Age ,data = test)

## Test Error ##
pred.final_mod <- predict(final.mod, newdata = test[,which(colnames(train) %in% c(sel.var.final, "BPSysAve"))], type = "response")
## Prediction error ##
pred.error.final <- mean((test$BPSysAve - pred.final_mod)^2)


## diagnostics checking for the final model ##

vif(final.mod)

## The hat values ###

h.final <- hatvalues(final.mod)
thresh <- 2 * (dim(model.matrix(final.mod))[2])/nrow(train)
w.final <- which(h.final > thresh)
w.final

## Standardized residuals ##

final_mod.stres <- rstudent(final.mod)

fitted.final <- predict(final.mod)
qqnorm(final_mod.stres)
plot(final_mod.stres ~ fitted.final, main= "Residual plot for the final model", type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

## Cook's distance##
D.final <- cooks.distance(final.mod)
D.final
which(D.final > qf(0.5, 5, 168-5))


## DFFITS ##
dfits.final <- dffits(final.mod)
which(abs(dfits.final) > 2*sqrt(5/168))

## DFBETAS ##
dfb.final <- dfbetas(final.mod)
which(abs(dfb.final[,2]) > 2/sqrt(168))


## ANOVA ##

anova(final.mod)


pred.error.AIC
pred.error.BIC
pred.error.lasso
pred.error.final

## diagnostics checking for the transformed model ##

vif(transform_mod)

## The hat values ###

h.transform <- hatvalues(transform_mod)
thresh <- 2 * (dim(model.matrix(transform_mod))[2])/nrow(train)
w.transform <- which(h.transform > thresh)
w.transform

## Standardized residuals ##

transform_mod.stres <- rstudent(transform_mod)

fitted.transform <- predict(transform_mod)
qqnorm(transform_mod.stres)
plot(transform_mod.stres ~ fitted.transform, main="Residual plot for the transformed model", type = "p", xlab = "Fitted Values", 
     ylab = "Standardized Residual", cex.lab = 1.5,
     col = "blue")

## Cook's distance##
D.transform <- cooks.distance(transform_mod)
D.transform
which(D.transform > qf(0.5, 5, 168-5))


## DFFITS ##
dfits.transform <- dffits(transform_mod)
which(abs(dfits.transform) > 2*sqrt(5/168))

## DFBETAS ##
dfb.transform <- dfbetas(transform_mod)
which(abs(dfb.transform[,2]) > 2/sqrt(168))


## ANOVA ##

anova(transform_mod)
