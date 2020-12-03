############################################################################
############################################################################
###########                                                      ###########
########### Script for plotting age prediction error on each     ###########
###########      dataset: balidation and test set and including  ###########
###########      statistical tests for age group comparisons     ###########
###########             Author: Diego Montiel Gonzalez           ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             d.montielgonzalez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

setwd("/PATH/Y-CpG/")


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


get.age.group <- function(df_test,pred){
  colnames(df_test_pred) <- c("Age", "AgePred")
  df_test_pred <- as.data.frame(df_test_pred)

  a <- df_test_pred[which(df_test_pred$Age <=20),]
  #a <- (a$Age - a$AgePred)
  a <- abs(a$Age - a$AgePred)
  b <- df_test_pred[which(df_test_pred$Age > 20 & df_test_pred$Age <= 40),]
  b <- abs(b$Age - b$AgePred)
  #b <- (b$Age - b$AgePred)
  c <- df_test_pred[which(df_test_pred$Age > 40 & df_test_pred$Age <= 60),]
  c <- abs(c$Age - c$AgePred)
  #c <- (c$Age - c$AgePred)
  d <- df_test_pred[which(df_test_pred$Age > 60),]
  d <- abs(d$Age - d$AgePred)
  #d <- (d$Age - d$AgePred)

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
  return(age_group)
}


library(Metrics)
library(ggplot2)

#################################################################
### Import data
#################################################################

df.train <- read.csv("data/normalized_Y-CpG/train.csv", row.names = "X")
df.validation <- read.csv("data/normalized_Y-CpG/validation.csv", row.names = "X")
df.test <- read.csv("data/normalized_Y-CpG/test.csv", row.names = "X")

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
### SVM regression
### Validation set
###############################################
svm.class <- svm(age ~ ., data = df.train.iqr, kernel = "radial", scale = TRUE, cost = 2, type = 'eps-regression')
pred <- as.numeric(predict(svm.class, df.val.iqr))

print(aad(df.validation.pheno$age - pred))
print(mae(df.validation.pheno$age, pred))
print(cor.test(df.validation.pheno$age, pred))
print(RMSE(df.validation.pheno$age, pred))
print(R2(df.validation.pheno$age, pred))

## Plotting Real Age vs Predicted Age
df.validation.pheno$pred <- pred
#pdf("validation_pred.pdf", width=10, height = 7)
qplot(df.validation.pheno$age, df.validation.pheno$pred, main="SVM radial validation") +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + ylim(13, 82) + xlim(13,82) +
  ylab("Predicted Age") + xlab("Real Age")   +
  geom_point(col = "bisque3", alpha = 0.9, size = 2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
#dev.off()

###############################################
### SVM regression
### Test set
###############################################
pred <- as.numeric(predict(svm.class, df.test.iqr))
mae(df.test.pheno$age, pred)
mdae(df.test.pheno$age, pred)
cor.test(df.test.pheno$age, pred)
RMSE(df.test.pheno$age, pred)
R2(df.test.pheno$age, pred)
df.test.pheno$pred <- pred

qplot(df.test.pheno$age, df.test.pheno$pred, main="SVM radial test") + geom_abline(intercept = 0, slope = 1, linetype = 2) + ylim(13, 82) + xlim(13,82) +
  ylab("Predicted Age") + xlab("Real Age")   +
  geom_point(col = "cornflowerblue", alpha = 0.9, size = 2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


##############################################
### Age groups and violin plots
##############################################
flag <- "validation"
#flag <- "test"

if(flag == "validation"){
  df_test_pred <- cbind(df.validation.pheno$age, df.validation.pheno$pred) #validation set
  values = c("#e3d8cd", "#d5c4b4", "#c0a78e", "#ab8a68")
}else if(flag == "test"){
  df_test_pred <- cbind(df.test.pheno$age, df.test.pheno$pred) # test set
  values <- c("#BFEFFF", "#9AC0CD", "#75A1D0")
}

age_group <- get.age.group(df_test_pred)

ggplot(age_group, aes(x=group, y=error, fill=group)) + labs(y="Prediction error in years") +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.65) +
  ggtitle("Errors dataset") +
  scale_fill_manual(values=values)  + # validation
  #scale_fill_manual(values=c("#BFEFFF", "#9AC0CD", "#75A1D0"))  + # testing
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_summary(fun.data=data_summary)


colnames(df_test_pred) <- c("Age", "AgePred")
df_test_pred <- as.data.frame(df_test_pred)

group.a <- df_test_pred[df_test_pred$Age <=20,]
mae(group.a$Age,group.a$AgePred)

group.b <- df_test_pred[df_test_pred$Age >20 & df_test_pred$Age <= 40,]
mae(group.b$Age,group.b$AgePred)

group.c <- df_test_pred[df_test_pred$Age >40 & df_test_pred$Age <= 60,]
mae(group.c$Age,group.c$AgePred)

group.d <- df_test_pred[df_test_pred$Age > 60,]
mae(group.d$Age,group.d$AgePred)

group.a$diff <- abs(group.a$Age - group.a$AgePred)
group.b$diff <- abs(group.b$Age - group.b$AgePred)
group.c$diff <- abs(group.c$Age - group.c$AgePred)
group.d$diff <- abs(group.d$Age - group.d$AgePred)

plot(density(group.a$diff))
plot(density(group.b$diff))
plot(density(group.c$diff))
plot(density(group.d$diff))
qqPlot(group.a$diff)
qqPlot(group.b$diff)
qqPlot(group.c$diff)
qqPlot(group.d$diff)

## Test for normality
## > 0.05 -> acceptes hypthesis of normally distributed
# validation set
p.vals <- c(shapiro.test(group.a$diff)$p.value, shapiro.test(group.b$diff)$p.value,
            shapiro.test(group.c$diff)$p.value, shapiro.test(group.d$diff)$p.value )
# testing set
p.vals <- c(shapiro.test(group.b$diff)$p.value, shapiro.test(group.c$diff)$p.value, shapiro.test(group.d$diff)$p.value )
p.adjust(p.vals, method="bonferroni")

########################################################################
## normally distributed
########################################################################
res.aov <- aov(error ~ group, data = age_group)
summary(res.aov)
plot(res.aov, 1)
plot(res.aov, 2)
TukeyHSD(res.aov)
pairwise.t.test(age_group$error, age_group$group, p.adjust.method = "BH")
leveneTest(error ~ group, data = age_group) # >0.05 -> variance statistically diff between groups
oneway.test(error ~ group, data = age_group)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

########################################################################
## non-parametric
########################################################################

kruskal.test(x = age_group$error, g = age_group$group)
pairwise.wilcox.test(age_group$error, age_group$group, p.adjust.method = "bonferroni")

