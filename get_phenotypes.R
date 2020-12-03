##############################################
### Get phenotypes using GEO query
##############################################

library(Biobase)
library(GEOquery)

getHeaderGeo <- function(accession)
{
  phenotype <- getGEO(accession, destdir=".")
  pheno <- phenotype[[1]]; rm(phenotype)
  pheno <- phenoData(pheno)
  pheno <- pData(pheno)
  filename <- paste0(accession,"_header.csv")
  write.csv(file = filename, pheno )
}

setwd("PATH/phenotypes/")

### "GSE128235"
accession <- "GSE128235"
#getHeaderGeo(accession)
gse.GSE128235 <- read.csv('GSE128235_header.csv')
gse.GSE128235 <- as.matrix(cbind(accession, unfactor(gse.GSE128235$geo_accession),  gse.GSE128235$age.ch1, unfactor(gse.GSE128235$Sex.ch1), unfactor(gse.GSE128235$tissue.ch1)))
colnames(gse.GSE128235) <- c("geo","accession","age","sex","tissue")

### "GSE100386"
accession <- "GSE100386"
#getHeaderGeo(accession)
gse.GSE100386 <- read.csv('GSE100386_header.csv')
gse.GSE100386 <- as.matrix(cbind(accession, unfactor(gse.GSE100386$geo_accession),  gse.GSE100386$age.ch1, unfactor(gse.GSE100386$Sex.ch1), unfactor(gse.GSE100386$source_name_ch1)))
colnames(gse.GSE100386) <- c("geo","accession","age","sex","tissue")

### "GSE125105"
accession <- "GSE125105"
#getHeaderGeo(accession)
gse.GSE125105 <- read.csv('GSE125105_header.csv')
gse.GSE125105 <- as.matrix(cbind(accession, unfactor(gse.GSE125105$geo_accession),  gse.GSE125105$age.ch1, unfactor(gse.GSE125105$Sex.ch1), unfactor(gse.GSE125105$tissue.ch1)))
colnames(gse.GSE125105) <- c("geo","accession","age","sex","tissue")

### "GSE61496"
accession <- "GSE61496"
#getHeaderGeo(accession)
gse.GSE61496 <- read.csv('GSE61496_header.csv')
gse.GSE61496 <- as.matrix(cbind(accession, unfactor(gse.GSE61496$geo_accession),  gse.GSE61496$age.ch1, gse.GSE61496$sex..1.m..2.f.ch1, unfactor(gse.GSE61496$tissue.ch1)))
colnames(gse.GSE61496) <- c("geo","accession","age","sex","tissue")
# sex, 1=m, 2=f: 2

### "GSE87571"
accession <- "GSE87571"
#getHeaderGeo(accession)
gse.GSE87571 <- read.csv('GSE87571_header.csv')
gse.GSE87571 <- as.matrix(cbind(accession, unfactor(gse.GSE87571$geo_accession),  gse.GSE87571$age.ch1, unfactor(gse.GSE87571$gender.ch1), unfactor(gse.GSE87571$tissue.ch1)))
colnames(gse.GSE87571) <- c("geo","accession","age","sex","tissue")

### "GSE115278"
accession <- "GSE115278"
#getHeaderGeo(accession)
gse.GSE115278 <- read.csv('GSE115278_header.csv')
colnames(gse.GSE115278)
gse.GSE115278 <- as.matrix(cbind(accession, unfactor(gse.GSE115278$geo_accession),  gse.GSE115278$age.ch1, unfactor(gse.GSE115278$Sex.ch1), unfactor(gse.GSE115278$source_name_ch1)))
colnames(gse.GSE115278) <- c("geo","accession","age","sex","tissue")

train.pheno <- as.data.frame(rbind(gse.GSE100386,gse.GSE125105,gse.GSE128235,gse.GSE61496,gse.GSE87571))
unique(as.vector(train.pheno$sex))
tissues.male <- c("Male","M","1")
#males = 1, M, Male
train.pheno <- train.pheno[which(train.pheno$sex %in% tissues.male),]
train.pheno$sex <- "M"

test.pheno <- as.data.frame(gse.GSE115278)
unique(test.pheno$sex)
test.pheno <- test.pheno[which(test.pheno$sex %in% tissues.male),]
test.pheno$sex <- "M"

train.bmiq <- fread("/media/disk1/diego/Epigenetics/Y-CpG/iDATs_Blood/df_bmiq_ycpgs_train.txt")
train.bmiq <- train.bmiq$Accession
validation.bmiq <- fread("/media/disk1/diego/Epigenetics/Y-CpG/iDATs_Blood/df_bmiq_ycpgs_test.txt")
validation.bmiq <- validation.bmiq$Accession
test.bmiq <- fread("/media/disk1/diego/Epigenetics/Y-CpG/iDATs_Blood/source/test/GSE115278/GSE115278_phenotypes.csv")
test.bmiq <- test.bmiq$Accession

df.train <- train.pheno[train.pheno$accession %in% train.bmiq, ]
df.validation <- train.pheno[train.pheno$accession %in% validation.bmiq, ]
df.test <- test.pheno[test.pheno$accession %in% test.bmiq, ]

write.csv("train.pheno.csv", x = df.train )
write.csv("validation.pheno.csv", x = df.validation )
write.csv("test.pheno.csv", x = df.test )

train.pheno$age <- unfactor(train.pheno$age)
test.pheno$age <- unfactor(test.pheno$age)

tmp.a <- cbind(df.train$age, "train")
tmp.b <- cbind(df.validation$age, "validation")
tmp.c <- cbind(df.test$age, "test")

df.age <- as.matrix(rbind(tmp.a, tmp.b, tmp.c))
colnames(df.age) <- c("age","dataset")
df.age <- as.data.frame(df.age)
color <- as.factor(df.age$dataset)
df.age$age <- as.numeric(df.age$age)

# "test"       "train"      "validation"
levels(color) <- c("cornflowerblue", "antiquewhite","bisque4")
color <- as.character(color)
table(color)

pdf(file = "/media/disk1/diego/Epigenetics/Y-CpG/iDATs_Blood/source/plots/age_dist.pdf", height = 7, width = 10)
h <- hist(unfactor(df.train$age), breaks = 72,   main="Blood 450k Males", xlab = "Age in years",xlim=c(15,87), xaxp=c(15,87,12), col= 'antiquewhite', border = 'antiquewhite')
h <- hist(unfactor(df.validation$age), breaks = 72,   main="Blood 450k Males", xlab = "Age in years",xlim=c(15,87), xaxp=c(15,87,12), col= 'bisque3', add = T, border = 'bisque3')
h <- hist(unfactor(df.test$age), breaks = 50,   main="Blood 450k Males", xlab = "Age in years",xlim=c(15,87), xaxp=c(15,87,12), col= alpha("cornflowerblue", 0.7), add = T, border = alpha("cornflowerblue", 0.01))
dev.off()

