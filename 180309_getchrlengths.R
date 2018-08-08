# @Author: Francis Blokzijl & Myrthe Jager
# @Date: 15 Mar 2017
# @Modified: 9 March 2018
# @Description: Get chromosome lengths of mm10 genome




#1 Define outdir
outdir = "~/surfdrive/Shared/ERCC1/Data/callable_genome/determine_callable/"

#2 Load mouse genome mm10
#2A Install & load BSgenome package
# biocLite("BSgenome")
library(BSgenome)
#2B select genome
ref_genome = "BSgenome.Mmusculus.UCSC.mm10" 
#2C Download & load genome
# biocLite(ref_genome)
library(ref_genome, character.only = TRUE)

#3 Autosomal seqlengths
#3A Calculate
autosomal_seqlengths = seqlengths(get(ref_genome))[1:19]
#3B Save
#write.table(autosomal_seqlengths, paste(outdir,"autosomal_chr_lengths_mm10.txt",sep=""), quote = F, sep = "\t", col.names=F)