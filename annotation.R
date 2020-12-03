############################################################################
############################################################################
###########                                                      ###########
###########       Get annotation for Y-CpGs using rgSet object   ###########
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


setwd("/PATH/Y-CpG/data/annotation/")

# load lbraries
library(minfi)
library(data.table)

features <- read.table("../feature_selection/IQR_features.txt")$V1
rgSet <- read.metharray.exp(getwd())
annotation  <- getAnnotation(rgSet) # get info from IDAT rgSet object
annotation <- annotation[features,]

igv_stochCpG <- data.frame(seqname = annotation$chr,
                           start = annotation$pos, end = annotation$pos+1,
                           name = rownames(annotation),
                           score = 1000)
write.table(x = 'track name=features description=Clone useScore=1', file = 'annotation_IGV.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = igv_stochCpG, file = 'annotation_IGV.bed', nThread = 4, sep = '\t', col.names = F, append = T)
