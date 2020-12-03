############################################################################
############################################################################
###########                                                      ###########
###########       Predict Age prediction using Y-CpG probes      ###########
###########                                                      ###########
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

#################################################################
### Import data
#################################################################

## Input Example load in a dataframe Normalized Y-CpGs

# Example with test set
dataset <- read.csv("data/normalized/normalized_Y-CpG/test.csv", row.names = "X")
df.test.pheno <- read.csv("data/test.pheno.csv")
dataset <- t(dataset) # transpose

# Your dataframe
# matrix: Samples as rows, CpGs as columns
dataset <- read.csv("PATH/YOUR_NORMALIZED_BETA_VALUES_YCPG.csv", row.names = "X")

## load full model
load("data/SVM_full.RData")
pred <- predict(svm.class, dataset)
pred
#get.metrics(df.test.pheno$age, pred) #in case of known age

## load reduced model
load("data/SVM_reduced.RData")
pred <- predict(svm.class, dataset)
pred
#get.metrics(df.test.pheno$age, pred)


