############################################################################
############################################################################
###########                                                      ###########
########### Script for data Quality control for IDATs (RAW 450K) ###########
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

setwd("PATH/Y-CpG/IDAT/")

## IDATs QC

# 0) Read IDATs in a rgSet object
rgSet <- read.metharray.exp(getwd(), extended = T)
qcinfo <- ENmix::QCinfo(rgSet)

# 1) Detect failed probes and samples with SQN w/rgSet
badsample <- unlist(lapply(X = strsplit(qcinfo$badsample, split = '[_]'), FUN = function(X) {X[1]})) # split the header names
length(qcinfo$badCpG)
write.csv(qcinfo$badsample, file="PATH/badsamples.txt", quote = F, row.names = F)
write.csv(qcinfo$badCpG, file="PATH/badCpGs.txt", quote = F, row.names = F)

# 2) Predict Sex on IDAT samples with Qantile Normalization
um.sqn <- preprocessQuantile(rgSet)
pdf(file = "sex_estimation_raw.pdf", width = 10, height = 7)
plotSex(um.sqn)
dev.off()

# 3) Plot density beta methylation values
pdf(file = "densityplot_raw.pdf", width = 10, height = 7)
densityPlot(rgSet)
dev.off()
