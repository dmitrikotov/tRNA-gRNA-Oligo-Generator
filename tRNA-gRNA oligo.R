#Script for designing oligos for cloning tRNA-gRNA arrays

library(Biostrings)

setwd("~/Documents/tRNA-gRNA Oligos")
guides <- read.csv("Jessica Guides 021220.csv")

Aar1_tRNA = "TATCACCTGCCCCC"
Aar1_rep  = "TATCACCTGCCCCA"
tRNA = "TGCACCAGCCGGGAATCG"
rep = "GTTTTAGAGCTAGAAATAGC"

oligos <- function(x){
x <- DNAString(x)
rep_guide = substr(x,9,20)
if (substr(as.character(x),1,4) == "GGTG") stop("Guide generates an AarI cut site")
tRNA_guide = reverseComplement(substr(x,1,12))
tRNA_oligo = paste(Aar1_tRNA,tRNA_guide,tRNA, sep = "")
rep_oligo = paste(Aar1_rep,rep_guide,rep, sep = "")
df <- data.frame(tRNA_oligo, rep_oligo)
return(df)
}

results <- as.data.frame(matrix(0, nrow(guides), 3))
colnames(results) <- c("gene","tRNA_oligo","rep_oligo")

for (i in 1:nrow(guides)) {
  results[i,2] <- as.character(oligos(guides[i,2])[1,1])
  results[i,3] <- as.character(oligos(guides[i,2])[1,2])
  results[i,1] <- as.character(guides[i,1])
}

write.csv(as.data.frame(results),file="tRNA-gRNA oligo.csv",quote=F)
