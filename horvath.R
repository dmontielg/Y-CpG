############################################################################
############################################################################
###########                                                      ###########
########### Horvath Clock comparison with Y-CpG model            ###########
###########             Author: Diego Montiel Gonzalez           ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             d.montielgonzalez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

# load libraries
library(wateRmelon)
library(data.table)
library(Metrics)
library(MLmetrics)

setwd("/PATH/Y-CpG")

## Input data
df.train.pheno <- read.csv("data/train.pheno.csv")
df.validation.pheno <- read.csv("data/validation.pheno.csv")
df.test.pheno <- read.csv("data/test.pheno.csv")

df.h.train <- read.csv("data/normalized_horvath/dataset_train_horvath.csv", row.names = "X")
df.h.val <- read.csv("data/normalized_horvath/dataset_validation_horvath.csv", row.names = "X")
df.h.test <- read.csv("data/normalized_horvath/dataset_test_horvath.csv", row.names = "X")

##############################
### HOrvarth Clock
##############################


data(coef) # load horvath coefficients

## training set
pred <- agep(df.h.train)[,1]
qplot(df.train.pheno$age, pred, main="Horvath clock in training set") + geom_abline(intercept = 0, slope = 1) + ylim(0,100) + xlim(0,100)
mae(pred, df.train.pheno$age)
cor(pred, df.train.pheno$age)

a <- cbind(df.train.pheno$age, pred)
a <- as.data.frame(a)
a["dataset"] <- "train"
colnames(a) <- c("age", "pred", "dataset")

## validation set
pred <- agep(df.h.val)[,1]
qplot(df.validation.pheno$age, pred, main="Horvath clock in validation set") + geom_abline(intercept = 0, slope = 1) + ylim(0,100) + xlim(0,100)
mae(pred, df.validation.pheno$age)
cor(pred, df.validation.pheno$age)

b <- cbind(df.validation.pheno$age, pred)
b <- as.data.frame(b)
b["dataset"] <- "validation"
colnames(b) <- c("age", "pred", "dataset")

## testing set
pred <- agep(df.h.test)[,1]
qplot(df.test.pheno$age, pred, main="Horvath clock in testing set") + geom_abline(intercept = 0, slope = 1) + ylim(0,100) + xlim(0,100)
mae(pred, df.test.pheno$age)
cor(pred, df.test.pheno$age)
RMSE(pred, df.test.pheno$age)
R2_Score(pred, df.test.pheno$age)

c <- cbind(df.test.pheno$age, pred)
c <- as.data.frame(c)
c["dataset"] <- "test"
colnames(c) <- c("age", "pred", "dataset")

horvath.dataframe <- rbind(a, b, c)

color <- horvath.dataframe$dataset
color <- as.factor(color)

# "test"       "train"      "validation"
levels(color) <- c("cornflowerblue", "antiquewhite","bisque4")
color <- as.character(color)
unique(color)
table(horvath.dataframe$dataset)


ggplot(horvath.dataframe, aes(x = age,y = pred, fill = dataset)) + ylim(c(13,82)) + xlim(13,82) +
  geom_point(alpha = 0.95, color = color) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

aad(horvath.dataframe$pred - horvath.dataframe$age)
mae(horvath.dataframe$pred, horvath.dataframe$age)
cor.test(horvath.dataframe$pred, horvath.dataframe$age)
cor.test(horvath.dataframe$pred, horvath.dataframe$age)$p.value
RMSE(horvath.dataframe$pred, horvath.dataframe$age)
R2_Score(horvath.dataframe$pred, horvath.dataframe$age)

df.h <- as.data.frame(cbind(df.train.pheno$age, t(df.h.train)))
df.h[1:10,1:10]
dim(df.h)
colnames(df.h)[1] <- "age"

#df_test_pred <- cbind(horvath.dataframe$age, horvath.dataframe$pred) #all dataset
#df_test_pred <- horvath.dataframe[horvath.dataframe$dataset == "validation",] #test dataset
df_test_pred <- horvath.dataframe[horvath.dataframe$dataset == "test",] #test dataset

colnames(df_test_pred) <- c("Age", "AgePred")
df_test_pred <- as.data.frame(df_test_pred)

a <- df_test_pred[which(df_test_pred$Age <=20),]
a <- abs(a$Age - a$AgePred)
b <- df_test_pred[which(df_test_pred$Age > 20 & df_test_pred$Age <= 40),]
b <- abs(b$Age - b$AgePred)
c <- df_test_pred[which(df_test_pred$Age > 40 & df_test_pred$Age <= 60),]
c <- abs(c$Age - c$AgePred)
d <- df_test_pred[which(df_test_pred$Age > 60),]
d <- abs(d$Age - d$AgePred)

tmp_a <- NULL
tmp_a$error <- as.matrix(as.data.frame(a))
tmp_a$group <- rep("<=20",length(a))
tmp_a <- as.data.frame(tmp_a)
colnames(tmp_a) <- c("error", "group")

tmp_b <- NULL
tmp_b$error <- as.matrix(as.data.frame(b))
tmp_b$group <- rep(">20-<=40",length(b))
tmp_b <- as.data.frame(tmp_b)
colnames(tmp_b) <- c("error", "group")

tmp_c <- NULL
tmp_c$error <- as.matrix(as.data.frame(c))
tmp_c$group <- rep(">40-<=60",length(c))
tmp_c <- as.data.frame(tmp_c)
colnames(tmp_c) <- c("error", "group")

tmp_d <- NULL
tmp_d$error <- as.matrix(as.data.frame(d))
tmp_d$group <- rep(">60",length(d))
tmp_d <- as.data.frame(tmp_d)
colnames(tmp_d) <- c("error", "group")

age_group <- NULL
age_group <- rbind(tmp_a,tmp_b,tmp_c,tmp_d)
table(age_group$group)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(age_group, aes(x=group, y=error, fill=group)) + labs(y="Prediction error in years") +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) +
  ggtitle("Errors test set") +
  scale_fill_manual(values=c("darkseagreen1", "darkseagreen2", "darkseagreen3", "darkseagreen4"))  +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_summary(fun.data=data_summary)

