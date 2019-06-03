# @Date: 27 march 2019
# @Author: Myrthe Jager
# @Modified: 
# @Description: Refit 96-type profiles with the new signatures
# Abbreviations: 



# ---- GET STARTED ----

#1 Load required packages
library(MutationalPatterns)
library("ggplot2")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(pracma)
library(reshape2)
library(NMF)
library("gridExtra")
library(grid)

#2 Define input and output directory
dir = ""
indir = paste(dir,"ERCC1/Data/",sep="")
outdir = paste(dir,"ERCC1/Results/SNV/newsigs/",sep="")

#3 Functions
# Fit to signatures with TSB
fit_to_signatures_s = function(mut_matrix, signatures)
{
  # make sure dimensions of input matrix are correct
  if (dim(mut_matrix)[1] != 192)
    stop( paste("Mutation count matrix input should have",
                "dimensions 96 X n samples") )
  
  if (dim(signatures)[1] != 192)
    stop("Signatures input should have dimensions 96 X n signatures")
  
  n_samples = dim(mut_matrix)[2]
  n_signatures = dim(signatures)[2]
  lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
  lsq_reconstructed = matrix(NA, nrow=192, ncol=n_samples)
  
  # Process each sample
  for (i in 1:ncol(mut_matrix))
  {
    y = mut_matrix[,i]
    lsq = lsqnonneg(signatures, y)
    lsq_contribution[,i] = lsq$x
    lsq_reconstructed[,i] = signatures %*% as.matrix(lsq$x) 
  }
  
  # Add row and col names
  sample_names = colnames(mut_matrix)
  signature_names = colnames(signatures)
  mut_type_names = rownames(signatures)
  
  colnames(lsq_contribution) = sample_names
  rownames(lsq_contribution) = signature_names
  
  colnames(lsq_reconstructed) = sample_names
  rownames(lsq_reconstructed) = mut_type_names
  
  res = list(lsq_contribution, lsq_reconstructed)
  names(res) = c("contribution", "reconstructed")
  
  return(res)
}

#4 Install and load mouse reference genome
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 Metadata
tissue = c(rep(c("Liver","Small intestine"),3),"Liver",rep(c("Liver","Small intestine"),2))
mouse = c(rep("Ercc1(-/D)",6),rep("WT",5))
type = c(rep(c("Ercc1(-/D) Liver","Ercc1(-/D) Small intestine"),3),"WT Liver",rep(c("WT Liver","WT Small intestine"),2))
mouse_order = c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")
mousetype = factor(type, levels = mouse_order)




# ---- GET VCFS ----

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




# ---- MUTATION SPECTRA ----

#1 Get 6 types
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

#2 Get 96 mutation types
mut_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

#3 Collapsed these profiles
collapsed_mutmatrix <- data.frame(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")

#4 Generate combined mutation numbers
collapsed_mutmatrix_combined <- data.frame(round(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1)/3),
                                           round(rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1)/3),
                                           round(rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1)/2),
                                           round(rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1)/3))
colnames(collapsed_mutmatrix_combined) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")

#5 Generate table with cosmic signatures
cancer_signatures_new = read.csv(paste(indir,"SNV/SBS_signatures.csv",sep=""))
cancer_signatures_new = cancer_signatures_new[order(cancer_signatures_new[,1]),]
#Tested and order of mutation types is similar to order in the old sigs
cancer_signatures_new = as.matrix(cancer_signatures_new[,3:67])

#5 Fit to sigs
fit_res_cancersigs = fit_to_signatures(collapsed_mutmatrix, cancer_signatures_new)
p1.1 <- plot_contribution(fit_res_cancersigs$contribution, coord_flip = T, signatures = cancer_signatures_new, mode = "absolute")
#ggsave(paste(outdir,"refit_all_newsigs.pdf",sep=""), p1.1, width = 12, height = 5)

#6 Define number of signatures
#6A
fit_res_collapsed_combined= fit_to_signatures(collapsed_mutmatrix_combined, cancer_signatures_new)
#6B Generate df
test.sigs <- data.frame(fit_res_collapsed_combined$contribution)
test.sigs$sum <- rowSums(test.sigs)
test.sigs$rank <- NA
test.sigs$sig <- NA
test.sigs$rnr <- NA
#6C Rank signatures, based on contribution
rank.sigs <- order(-test.sigs$sum)
for(i in rank.sigs) {
  sig <- rank.sigs[i]
  test.sigs[sig,]$sig <- substring(rownames(test.sigs[sig,]), 4)
  
  test.sigs[sig,]$rnr <- as.numeric(sig)
  
  rank <- i
  test.sigs[sig,]$rank <- rank
  remove(sig,rank)
}
test.sigs <- test.sigs[order(-test.sigs$sum),]
#6D Generate new empty df
df.sigs.test <- data.frame(n.sigs = NA,
                           wt.liver = NA,
                           mut.liver = NA,
                           wt.si = NA,
                           mut.si = NA,
                           mean = NA)
#6E Calculate how much of the original profile is reconstructed by the refitted data
# First 2 sigs, then 3, etc.
for(i in 2:nrow(test.sigs)) {
  sigs.to.test <- test.sigs[1:i,]$rnr
  cancer.sigs <- cancer_signatures_new[,sigs.to.test]
  reconstructed.df <- fit_to_signatures(collapsed_mutmatrix,cancer.sigs)$reconstructed
  wt.liver <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix)[1,1]
  mut.liver <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix)[2,2]
  wt.si <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix)[3,3]
  mut.si <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix)[4,4]
  mean.all <- mean(wt.liver,mut.liver,wt.si,mut.si)
  df.sigs.test <- rbind(df.sigs.test,c(i,wt.liver,mut.liver,wt.si,mut.si, mean.all))
  remove(sigs.to.test,cancer.sigs,wt.liver,mut.liver,wt.si,mut.si, mean.all)
}
df.sigs.test <- df.sigs.test[-1,]
#6F Plot
df.sigs.test.forplot <- melt(df.sigs.test[,-6], id.vars = "n.sigs")
p1.2 <- ggplot(df.sigs.test.forplot, aes(x=n.sigs,y=value,fill= variable)) +
  geom_bar(stat = "identity", position = position_dodge())+ 
  ggtitle("Collapsed & adjusted for n mouse") +
  geom_hline(yintercept = 0.90) +
  scale_x_continuous("Signatures", labels = as.character(2:65), breaks = c(2:65)) +
  coord_cartesian(ylim=c(0.8, 1)) +
  theme_bw()
#ggsave(paste(outdir,"nsigs_collapsed_newsigs.pdf",sep=""),plot = p1.2, width = 12, height = 4)

#7 Choose threshold. When nothing changes; 25 signatures
impt.sigs_new <- sort(test.sigs[1:25,]$rnr)

#7 Fit to important sigs
fit_res_cancersigs_impt = fit_to_signatures(collapsed_mutmatrix, cancer_signatures_new[,impt.sigs_new])
p1.3 <- plot_contribution(fit_res_cancersigs_impt$contribution, coord_flip = T, signatures = cancer_signatures_new[,impt.sigs_new], mode = "absolute")
#ggsave(paste(outdir,"refit_impt_newsigs.pdf",sep=""), p1.3, width = 12, height = 5)




# ---- TRANSCRIPTIONAL STRAND BIAS ----

#1 Get knowngenes table from UCSC for genome
genes_mm10 = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

#2 Make mutation count matrix with transcriptional information
mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_mm10)


#3 Collapsed these profiles
collapsed_mutmatrix_s <- data.frame(rowSums(mut_mat_s[,c(7,8,10)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_mat_s[,c(1,3,5)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_mat_s[,c(9,11)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_mat_s[,c(2,4,6)], na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix_s) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")

#4 Generate combined mutation numbers
collapsed_mutmatrix_combined_s <- data.frame(round(rowSums(mut_mat_s[,c(7,8,10)], na.rm = FALSE, dims = 1)/3),
                                             round(rowSums(mut_mat_s[,c(1,3,5)], na.rm = FALSE, dims = 1)/3),
                                             round(rowSums(mut_mat_s[,c(9,11)], na.rm = FALSE, dims = 1)/2),
                                             round(rowSums(mut_mat_s[,c(2,4,6)], na.rm = FALSE, dims = 1)/3))
colnames(collapsed_mutmatrix_combined_s) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")

#5 Generate table with cosmic signatures
cancer_signatures_new_s = read.csv(paste(indir,"SNV/TSB_signatures.csv",sep=""))
cancer_signatures_new_s = cancer_signatures_new_s[order(cancer_signatures_new_s[,2],
                                                        cancer_signatures_new_s[,3]),]
#Tested and order of mutation types is similar to order in the old sigs
cancer_signatures_new_s = as.matrix(cancer_signatures_new_s[,4:68])

#6 Fit to sigs
fit_res_cancersigs_s = fit_to_signatures_s(collapsed_mutmatrix_s, cancer_signatures_new_s)
p2.1 <- plot_contribution(fit_res_cancersigs_s$contribution, coord_flip = T, signatures = cancer_signatures_new_s, mode = "absolute")
#ggsave(paste(outdir,"TSB_refit_all_newsigs.pdf",sep=""), p2.1, width = 12, height = 5)

#7 Define number of signatures
#7A
fit_res_collapsed_combined_s= fit_to_signatures_s(collapsed_mutmatrix_combined_s, cancer_signatures_new_s)
#7B Generate df
test.sigs_s <- data.frame(fit_res_collapsed_combined_s$contribution)
test.sigs_s$sum <- rowSums(test.sigs_s)
test.sigs_s$rank <- NA
test.sigs_s$sig <- NA
test.sigs_s$rnr <- NA
#7C Rank signatures, based on contribution
rank.sigs_s <- order(-test.sigs_s$sum)
for(i in rank.sigs_s) {
  sig_s <- rank.sigs_s[i]
  test.sigs_s[sig_s,]$sig <- substring(rownames(test.sigs_s[sig_s,]), 4)
  
  test.sigs_s[sig_s,]$rnr <- as.numeric(sig_s)
  
  rank_s <- i
  test.sigs_s[sig_s,]$rank <- rank_s
  remove(sig_s,rank_s)
}
test.sigs_s <- test.sigs_s[order(-test.sigs_s$sum),]
#7D Generate new empty df
df.sigs.test_s <- data.frame(n.sigs = NA,
                           wt.liver = NA,
                           mut.liver = NA,
                           wt.si = NA,
                           mut.si = NA,
                           mean = NA)
#7E Calculate how much of the original profile is reconstructed by the refitted data
# First 2 sigs, then 3, etc.
for(i in 2:nrow(test.sigs_s)) {
  sigs.to.test <- test.sigs_s[1:i,]$rnr
  cancer.sigs <- cancer_signatures_new_s[,sigs.to.test]
  reconstructed.df <- fit_to_signatures_s(collapsed_mutmatrix_s,cancer.sigs)$reconstructed
  wt.liver <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix_s)[1,1]
  mut.liver <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix_s)[2,2]
  wt.si <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix_s)[3,3]
  mut.si <- cos_sim_matrix(reconstructed.df,collapsed_mutmatrix_s)[4,4]
  mean.all <- mean(wt.liver,mut.liver,wt.si,mut.si)
  df.sigs.test_s <- rbind(df.sigs.test_s,c(i,wt.liver,mut.liver,wt.si,mut.si, mean.all))
  remove(sigs.to.test,cancer.sigs,wt.liver,mut.liver,wt.si,mut.si, mean.all)
}
df.sigs.test_s <- df.sigs.test_s[-1,]
#7F Plot
df.sigs.test.forplot_s <- melt(df.sigs.test_s[,-6], id.vars = "n.sigs")
p2.2 <- ggplot(df.sigs.test.forplot_s, aes(x=n.sigs,y=value,fill= variable)) +
  geom_bar(stat = "identity", position = position_dodge())+ 
  ggtitle("Collapsed & adjusted for n mouse") +
  geom_hline(yintercept = 0.90) +
  scale_x_continuous("Signatures", labels = as.character(2:65), breaks = c(2:65)) +
  coord_cartesian(ylim=c(0.8, 1)) +
  theme_bw()
#ggsave(paste(outdir,"TSB_nsigs_collapsed_newsigs.pdf",sep=""),plot = p2.2, width = 12, height = 4)

#8 Choose threshold. When nothing changes; 25 signatures
impt.sigs_new_s <- sort(test.sigs_s[1:23,]$rnr)

#7 Fit to important sigs
fit_res_cancersigs_impt_s = fit_to_signatures_s(collapsed_mutmatrix_s, cancer_signatures_new_s[,impt.sigs_new_s])
p2.3 <- plot_contribution(fit_res_cancersigs_impt_s$contribution, coord_flip = T, signatures = cancer_signatures_new_s[,impt.sigs_new_s], mode = "absolute")
#ggsave(paste(outdir,"TSB_refit_impt_newsigs.pdf",sep=""), p2.3, width = 12, height = 5)
