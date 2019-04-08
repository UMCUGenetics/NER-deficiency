# @Date: 24 october 2017
# @Modified: 08 august 2018
# @Author: Francis Blokzijl and Myrthe Jager
# @Description: Determine expression of genes with mutations




# ---- GET STARTED ----

#1 Load required packages
library(BSgenome)
library(MutationalPatterns)
library(org.Mm.eg.db)
library(plyr)
library(ggplot2)

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/SNV/"

#3 Functions
## Function 1: Get coding mutations ##
get_coding_muts <- function(vcflist) {
  overlap_total <- data.frame()
  for(i in 1:length(vcflist) )  {
    overlap <-data.frame()
    overlap = findOverlaps(vcflist[[i]], genes_mm10)
    overlap <- as.data.frame(overlap)
    colnames(overlap) <- c("mutations","genes")
    overlap$gene_id <- genes_mm10[overlap$genes]$gene_id
    overlap$sample <- names(vcflist[i])
    ranges.df <- as.data.frame(ranges(vcflist[[i]]))
    ranges.df <- ranges.df[overlap$mutations,c(1,4)]
    colnames(ranges.df) <- c("position","mutationname")
    overlap <- cbind(overlap,ranges.df)
    overlap_total <- rbind(overlap_total, overlap)
  }
  return(overlap_total)
}

## Function 2: get ENSEMBL id ##
get_ensembl_id <- function(genes) {
  # voeg symbol toe aan je genes object
  genes$symbol <- mapIds(org.Mm.eg.db,
                         keys=genes$gene_id,
                         column="SYMBOL",
                         keytype="ENTREZID",
                         multiVals="first")
  # en ensembl gene id
  genes$ensembl <- mapIds(org.Mm.eg.db,
                          keys=genes$gene_id,
                          column="ENSEMBL",
                          keytype="ENTREZID",
                          multiVals="first")
  genes <- genes[-which(is.na(genes$ensembl)),]
  return(genes)
}

## Function 3: Calculate median RPKM
get_median_RPKM <- function(df) {
  df_new <- as.data.frame(matrix(ncol = 2, nrow = 0))
  for (i in 1:nrow(df)) {
    median <- median(df[i,])
    average <- sum(df[i,])/length(df[i,])
    df_new <- rbind(df_new,c(median,average))
  }
  colnames(df_new) <- c("median","average")
  df <- cbind(df,df_new)
  return(df)
}

#4 Install and load mouse reference genome
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 Get mm10 txdb genes
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
genes_mm10 = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

#6 Metadata
tissue = c(rep(c("Liver","Small intestine"),3),"Liver",rep(c("Liver","Small intestine"),2))
mouse = c(rep("Ercc1(-/D)",6),rep("WT",5))
type = c(rep(c("Ercc1(-/D) Liver","Ercc1(-/D) Small intestine"),3),"WT Liver",rep(c("WT Liver","WT Small intestine"),2))
mouse_order = c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")
mousetype = factor(type, levels = mouse_order)




# ------ GET VCFS ------

#1 Get VCFS
vcf_files = list.files(paste(indir,"SNV/VCF/",sep=""), full.names = T)
sample_names = c("Ercc1(-/D)1 Liver", "Ercc1(-/D)1 Small intestine", 
                 "Ercc1(-/D)2 Liver", "Ercc1(-/D)2 Small intestine", 
                 "Ercc1(-/D)3 Liver", "Ercc1(-/D)3 Small intestine",  
                 "WT1 Liver", 
                 "WT2 Liver","WT2 Small intestine",
                 "WT3 Liver", "WT3 Small intestine")
vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)

#2 Check if chromosome names are uniform
all(seqlevels(vcfs[[1]]) %in% seqlevels(get(ref_genome)))

#3 only select autosomal chromosomes
auto = extractSeqlevelsByGroup(species="Mus_musculus", style="UCSC", group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))




# ---- GET MUTATIONS WITHIN GENES ----

#1 Generate separate lists per mousetype & tissue
vcfs_mutliver <- vcfs[c(1,3,5)]
vcfs_mutsi <- vcfs[c(2,4,6)]
vcfs_wtliver <- vcfs[c(7,8,10)]
vcfs_wtsi <- vcfs[c(9,11)]

#2 Get coding mutations
coding_mutliver <- get_coding_muts(vcfs_mutliver)
coding_mutsi <- get_coding_muts(vcfs_mutsi)
coding_wtliver <- get_coding_muts(vcfs_wtliver)
coding_wtsi <- get_coding_muts(vcfs_wtsi)

#3 Add gene IDs
coding_mutliver<- get_ensembl_id(coding_mutliver)
coding_mutsi <- get_ensembl_id(coding_mutsi)
coding_wtliver <- get_ensembl_id(coding_wtliver)
coding_wtsi <- get_ensembl_id(coding_wtsi)




# ---- Expression of mutated genes ----

#1 Load RPKM file
RPKM = read.delim(paste(indir,"RNAseq/170117_NS500413_0253_AHC5NTBGX2_MOUSE_readCounts_RPKM.txt",sep=""), header = T)

#2 Get RPKMs of genes per mousetype and tissue
RPKM_mutliver <- as.matrix(sapply(RPKM[,c(2,4,6)], as.numeric))
RPKM_mutsi <- as.matrix(sapply(RPKM[,c(3,5,7)], as.numeric))
RPKM_wtliver <- as.matrix(sapply(RPKM[,c(8,9,12)], as.numeric))
RPKM_wtsi <- as.matrix(sapply(RPKM[,c(10,11,13,14)], as.numeric))

#3 Calculate median RPKM
RPKM_mutliver <- get_median_RPKM(RPKM_mutliver)
RPKM_mutsi <- get_median_RPKM(RPKM_mutsi)
RPKM_wtliver <- get_median_RPKM(RPKM_wtliver)
RPKM_wtsi <- get_median_RPKM(RPKM_wtsi)

#4 Add ensembl id
RPKM_mutliver$ensembl <- as.character(RPKM[,1])
RPKM_mutsi$ensembl <- as.character(RPKM[,1])
RPKM_wtliver$ensembl <- as.character(RPKM[,1])
RPKM_wtsi$ensembl <- as.character(RPKM[,1])

#5 combine tables
mutliver <- merge(coding_mutliver,RPKM_mutliver)
mutsi <- merge(coding_mutsi,RPKM_mutsi)
wtliver <- merge(coding_wtliver,RPKM_wtliver)
wtsi <- merge(coding_wtsi,RPKM_wtsi)
df_mutliver <-  ddply(mutliver, .(mutationname), function(x) x[which.max(x$median),])
df_mutsi <-  ddply(mutsi, .(mutationname), function(x) x[which.max(x$median),])
df_wtliver <-  ddply(wtliver, .(mutationname), function(x) x[which.max(x$median),])
df_wtsi <-  ddply(wtsi, .(mutationname), function(x) x[which.max(x$median),])

#6 Generate new dfs with essential information
df_mutliver_final <- data.frame(df_mutliver$sample,df_mutliver$median,df_mutliver$average,"mutliver")
colnames(df_mutliver_final) <- c("animal","median","average","animaltype")
df_mutsi_final <- data.frame(df_mutsi$sample,df_mutsi$median,df_mutsi$average,"mutsi")
colnames(df_mutsi_final) <- c("animal","median","average","animaltype")
df_wtliver_final <-data.frame(df_wtliver$sample,df_wtliver$median,df_wtliver$average,"wtliver")
colnames(df_wtliver_final) <- c("animal","median","average","animaltype")
df_wtsi_final <- data.frame(df_wtsi$sample,df_wtsi$median,df_wtsi$average,"wtsi")
colnames(df_wtsi_final) <- c("animal","median","average","animaltype")

#7 Combine into a single table
df.final <- rbind(df_mutliver_final,df_mutsi_final,df_wtliver_final,df_wtsi_final)
df.final$animaltype = factor(df.final$animaltype, levels = c("wtliver","mutliver","wtsi","mutsi" ))
df.final$animal = factor(df.final$animal, levels = c ("WT1 Liver","WT2 Liver","WT3 Liver",
                             "Ercc1(-/D)1 Liver","Ercc1(-/D)2 Liver","Ercc1(-/D)3 Liver",
                             "WT2 Small intestine","WT3 Small intestine",
                             "Ercc1(-/D)1 Small intestine","Ercc1(-/D)2 Small intestine","Ercc1(-/D)3 Small intestine"
                             ))

#8 Plot RPKM
plot <- ggplot(df.final, aes(x = animaltype, y = average)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()

#9 Save plot
ggsave(paste(outdir,"averageexpression_geneswithmut_pertype.pdf",sep=""), plot = plot, width = 5, height = 5)

#10 T test
p1 <- t.test(df.final[which(df.final$animaltype == "mutliver"),]$average,df.final[which(df.final$animaltype == "wtliver"),]$average)$p.value
p2 <- t.test(df.final[which(df.final$animaltype == "mutsi"),]$average,df.final[which(df.final$animaltype == "wtsi"),]$average)$p.value
p.adjust(c(p1,p2), method = "fdr") # 0.5007507 0.1072669