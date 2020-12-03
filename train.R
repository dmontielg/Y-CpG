############################################################################
############################################################################
###########                                                      ###########
########### Script training model for Age prediction with SVM    ###########
###########             with and without feature selection       ###########
###########             Author: Diego Montiel Gonzalez           ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             d.montielgonzalez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

setwd("/PATH/Y-CpG")

get.metrics <- function(real, pred)
{
  print(paste0("MAD: ", aad(real - pred)))
  print(paste0("MAE: ",mae(real, pred)))
  print(cor.test(real, pred))
  print(paste0("RMSE: ", RMSE(real, pred)))
  print(paste0("R2: ",R2(real, pred)))
}


library(car)
library(Metrics)
library(varhandle)
library(lsr)
library(data.table)
library(e1071)
library(ggplot2)
library(caret)
library(matrixStats)
library(leaps)

#################################################################
### Import data
#################################################################

df.train <- read.csv("data/normalized/normalized_Y-CpG/train.csv", row.names = "X")
df.validation <- read.csv("data/normalized/normalized_Y-CpG/validation.csv", row.names = "X")
df.test <- read.csv("data/normalized/normalized_Y-CpG/test.csv", row.names = "X")

df.train.pheno <- read.csv("data/train.pheno.csv")
df.validation.pheno <- read.csv("data/validation.pheno.csv")
df.test.pheno <- read.csv("data/test.pheno.csv")

y.cpg <- as.vector(read.csv("data/Y_cpgs_300.txt", header = F)$V1) # list filtered by cross-reactives probes and SNPs
y.cpg.intersect <- intersect(y.cpg, rownames(df.train))
df.test <- df.test[y.cpg.intersect,]

#################################################################
### filtering by IQR
#################################################################
iqr.threshold <- 0.1
iqrs <- rowIQRs(as.matrix(df.train))
hist(iqrs, breaks = 'fd')
abline(v = iqr.threshold)
IQR.cpgs <- rownames(df.train)[which(iqrs >= iqr.threshold)]
length(IQR.cpgs)

df.train.iqr <- as.data.frame(t(df.train[IQR.cpgs,]))
df.train.iqr$age <- df.train.pheno$age
df.val.iqr <- as.data.frame(t(df.validation))
df.test.iqr <- as.data.frame(t(df.test[IQR.cpgs,]))

###############################################
### SVM regression with 75 CpGs
###############################################
svm.class <- svm(age ~ ., data = df.train.iqr, kernel = "radial", scale = TRUE, cost = 2, type = 'eps-regression')
save(svm.class, file = "data/SVM_full.RData")
### Validation set
val.pred <- predict(svm.class, df.val.iqr)
get.metrics(df.validation.pheno$age, val.pred)

###############################################
### SVM regression
### Testing set
###############################################
test.pred <- predict(svm.class, df.test.iqr)
get.metrics(df.test.pheno$age, test.pred)

##########################################
### Stepwise-Forward Feature Selection ###
##########################################
n = 75
subset.fwd <- regsubsets(age ~., data = df.train.iqr, nvmax = n, really.big = T, method = "forward")
coef(subset.fwd, n)
reg.summary <- summary(subset.fwd)
names(reg.summary)

plot(reg.summary$adjr2, xlab = "Number of features", ylab = "Adj. Rsq", type = "l")
r2 <- which.max(reg.summary$adjr2)
points(r2, reg.summary$adjr2[r2], col="red", cex=2, pch=20)

plot(reg.summary$cp, xlab = "Number of features", ylab = "Cp", type = "l")
aic_cp <- which.min(reg.summary$cp)
points(aic_cp, reg.summary$cp[aic_cp], col="red", cex=2, pch=20)

plot(reg.summary$bic, xlab = "Number of features", ylab = "BIC", type = "l")
reg.summary$bic
bay_ic <- which.min(reg.summary$bic)
points(bay_ic, reg.summary$bic[bay_ic], col="red", cex=2, pch=20)
bay_ic
coef(subset.fwd, bay_ic)
features_bic <- coef(subset.fwd, bay_ic)
#write.table(coef(subset.fwd, bay_ic), file = "../update/bic_subset_fwd.csv", quote = F)
feature_selection <- names(features_bic[2:length(features_bic)])
df.train.iqr.reduced <- cbind(df.train.iqr[feature_selection], df.train.iqr["age"] )
###############################################
### SVM regression
#### with Feature selection
###############################################
svm.class <- svm(age ~., data = df.train.iqr.reduced, kernel = "radial", scale = TRUE, cost = 2, type = 'eps-regression')
save(svm.class, file = "data/SVM_reduced.RData")

### Validation set
val.pred <- predict(svm.class, df.val.iqr)
get.metrics(df.validation.pheno$age, val.pred)
### Testing set with Feature selection
test.pred <- predict(svm.class, df.test.iqr)
get.metrics(df.test.pheno$age, test.pred)
