############################################################################
############################################################################
###########                                                      ###########
########### Script for data Normalization IDATs BMIQ + ENmix     ###########
###########             Author: Diego Montiel Gonzalez           ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             d.montielgonzalez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################

## Load libraries ##

library(minfi)
library(data.table)
library(wateRmelon)
library(ENmix)
library(RPMM)

## Path for IDATs after QC
setwd("PATH/Y-CpG/IDAT/")

#####################################
### Normalization with BMIQ + ENmix
#####################################
# train + validation IDATs
probes2remove <- read.table("../data/qc_train/probes2remove.txt",header = F)
probes2remove <- probes2remove$V1
probes2remove <- (probes2remove[!duplicated(probes2remove)])

## Tested with 40 cores ~2 hours
gc()
rgSet <- read.metharray.exp(getwd())
um.q0 <- preprocessENmix(rgSet, bgParaEst = 'oob', dyeCorr = 'RELIC',
                         nCores=40, exCpG = probes2remove)
qc();
# Save an object to a file
saveRDS(um.q0, file = "../data/um.q0_preprocessENMIx.rds")
# Restore the object
um.q0 <- readRDS(file = "../data/um.q0_preprocessENMIx.rds")
gc();
um.qi <- norm.quantile(mdat = um.q0, method='quantile1')
saveRDS(um.qi, file = "../data/um.qi_preprocessENMIx.rds")
rm(um.q0) ; gc()
gc()
um.qi <- readRDS(file = "../data/um.qi_preprocessENMIx.rds")

df.beta.bmiq <- ENmix::bmiq.mc(mdat = um.qi, nCores=40)
samples <- unlist(lapply(X = strsplit(colnames(df.beta.bmiq), split = '[_]'), FUN = function(X) {X[1]}))
fwrite(as.data.frame(nocombat.beta.bmiq), file = "../data/beta_values_bmiq.txt", sep="\t", nThread = 40, row.names = T)
