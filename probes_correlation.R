############################################################################
############################################################################
###########                                                      ###########
###########       Get correlations of Y-CpGs methylation levels  ###########
###########       with age                                       ###########
###########                    from minfi library                ###########
###########             Author: Diego Montiel Gonzalez           ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             d.montielgonzalez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################


setwd("/media/disk1/diego/git/Y-CpG")

# load libraries

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(data.table)
library(varhandle)

process.beta.fread <- function(beta)
{
  cpgs <- as.vector(beta$V1)
  beta <- beta[,-1]
  beta <- as.matrix(beta)
  rownames(beta) <- cpgs
  return(beta)
}

############################################
### Correlation Y-CpGs
############################################
df.train <- read.csv("data/normalized_Y-CpG/train.csv", row.names = "X")
df.validation <- read.csv("data/normalized_Y-CpG/validation.csv", row.names = "X")
df.test <- read.csv("data/normalized_Y-CpG/test.csv", row.names = "X")

df.train.pheno <- read.csv("data/train.pheno.csv")
df.validation.pheno <- read.csv("data/validation.pheno.csv")
df.test.pheno <- read.csv("data/test.pheno.csv")

y.cpg <- as.vector(read.csv("data/Y_cpgs_300.txt", header = F)$V1) # list filtered by cross-reactives probes and SNPs
y.cpg.intersect <- intersect(y.cpg, rownames(df.train))
df.test <- df.test[y.cpg.intersect,]

df.train.ycpg <- df.train[y.cpg.intersect,]
df.test.ycpg <- df.test[rownames(df.train.ycpg),]
df.train.ycpg <- t(df.train.ycpg)
df.test.ycpg <- t(df.test.ycpg)

df.pheno <- rbind(df.train.pheno, df.validation.pheno, df.test.pheno)
df.all <- rbind(df.train.ycpg, df.test.ycpg)
rownames(df.pheno) <- df.pheno$accession
df.pheno <- df.pheno[rownames(df.all),]

df <- NULL
age.all <- df.pheno$age

for(i in 1:ncol(df.all))
{
  cor.pvalue <- cor.test(as.numeric(df.all[,i]), age.all, method="spearman", exact=F)
  tmp <- data.frame(cbind(colnames(df.all)[i], as.numeric(cor.pvalue$estimate), as.numeric(cor.pvalue$p.value)))
  rownames(tmp) <- i
  df <- rbind(df,tmp)
}

df <- as.data.frame(df)
colnames(df) <- c("cpg","estimate","p-value")
rownames(df) <- df$cpg
df$cpg <- NULL
df$adj.pvalue.bonferroni <- p.adjust(as.vector(df$`p-value`), method = "bonferroni", n = nrow(df))
# write.csv(df, file="train_test_all-Ycpgs_correlation.txt")

### Here needs to be load the entire autosomal beta values for comparisons autosomal and Y-chrosomoe probes
#df.bmiq.autosomal <- as.matrix(df.train[setdiff(rownames(df.train), y.cpg),])
#dim(df.bmiq.autosomal)
#iqrs1 <- rowIQRs(df.bmiq.autosomal)
#iqrs2 <- colIQRs(df.train.ycpg)

#plot(density(iqrs1), xlab = "IQRs"); lines(density(iqrs2), col = "red")
#plot(ecdf(iqrs1), xlab = "IQRs", main="IQR cumulative distribution"); lines(ecdf(iqrs2), col = "red")
#legend(x = "right", legend = c("Autosomal", "Y-Chr"), col = c("black", "red"), pch = 19, bty = "n", cex = 1.5)
#text(x = 0.2, y = 0.6, labels = "p-value = 2.641e-06\n(Mann-Whitney U-test)")
# independent 2-group Mann-Whitney U Test
#wilcox.test(x = iqrs2, y = iqrs1, alternative = "greater") # where y and x are numeric

#mean(iqrs1) #[1] 0.04892996
#mean(iqrs2) #[1] 0.06916646

############################################
### Annotation and correlation: IGV
############################################
setwd("/PATH/Y-CpG/data/annotation/")


data(Locations)
y_cpgs <- as.vector(read.csv("../Y_cpgs_300.txt", header = F)$V1)

Locations = Locations[Locations$chr == "chrY",]
Locations = Locations[order(Locations$pos),]
Locations <- Locations[rownames(df),]

score <- scales::rescale(unfactor(df$estimate), to=c(0,1))*1000
#score = abs(rnorm(nrow(Locations), 500)) # Replace for correlation. IGV prefers scores between 0 and 1000; if not, change viewLimits in the IGV header
coverage_track <- data.frame(seqname = Locations$chr,
                             start = as.integer(Locations$pos), end = as.integer(Locations$pos),
                             name = paste('chrY', paste(Locations$pos, Locations$pos, sep = "-"), sep = ':'),
                             score = score)
# Write header
#write.table(x = 'track name=y_track description=Y_track useScore=1 graphType=line viewLimits=0:1000', file = 'y_track.bed',
#            quote = F, sep = '\t', row.names = F, col.names = F)

# Append data
#fwrite(x = coverage_track, file = 'y_track.bed', nThread = 4, sep = '\t', col.names = F, append = T)
#View(cbind(as.numeric(df$estimate), (as.numeric(df$estimate) * 1000)))
