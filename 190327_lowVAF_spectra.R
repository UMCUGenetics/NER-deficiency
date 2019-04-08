# @Date: 27 march 2019
# @Author: Myrthe Jager
# @Modified: 
# @Description: Plot mutational spectra of variants with VAF < 0.3
# Abbreviations: 




# ---- GET STARTED ----

#1 Load required packages
library(MutationalPatterns)
library("ggplot2")

#2 Define input and output directory
indir = "C:/Users/myrth/surfdrive/Shared/ERCC1/Data/" 
outdir = "C:/Users/myrth//surfdrive/Shared/ERCC1/Results/SNV/VAFbelow03/"

#3 Functions

#4 Install and load mouse reference genome
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 Metadata
tissue = rep(c(rep(c("Liver","Small intestine"),3),"Liver",rep(c("Liver","Small intestine"),2)),2)
mouse = rep(c(rep("Ercc1(-/D)",6),rep("WT",5)),2)
type = rep(c(rep(c("Ercc1(-/D) Liver","Ercc1(-/D) Small intestine"),3),"WT Liver",rep(c("WT Liver","WT Small intestine"),2)),2)
muttype = c(rep("SNV",11),rep("low VAF",11))
plottype = c(c(rep(c("Ercc1(-/D) Liver","Ercc1(-/D) Small intestine"),3),"WT Liver",rep(c("WT Liver","WT Small intestine"),2)),rep("low VAF",11))
type = paste(type,muttype)
mouse_order = c("WT Liver SNV","Ercc1(-/D) Liver SNV","WT Small intestine SNV","Ercc1(-/D) Small intestine SNV",
                "WT Liver low VAF","Ercc1(-/D) Liver low VAF","WT Small intestine low VAF","Ercc1(-/D) Small intestine low VAF")
mousetype = factor(type, levels = mouse_order)




# ---- GET VCFS ----

#1 Get VCFS
vcf_files_lowVAF = list.files(paste(indir,"SNV/SNVbelowVAF03/",sep=""), full.names = T)
sample_names_lowVAF = c("Ercc1(-/D)1 Liver_lowVAF", "Ercc1(-/D)1 Small intestine_lowVAF", 
                 "Ercc1(-/D)2 Liver_lowVAF", "Ercc1(-/D)2 Small intestine_lowVAF", 
                 "Ercc1(-/D)3 Liver_lowVAF", "Ercc1(-/D)3 Small intestine_lowVAF",  
                 "WT1 Liver_lowVAF", 
                 "WT2 Liver_lowVAF","WT2 Small intestine_lowVAF",
                 "WT3 Liver_lowVAF", "WT3 Small intestine_lowVAF")
vcf_files = list.files(paste(indir,"SNV/VCF/",sep=""), full.names = T)
sample_names = c("Ercc1(-/D)1 Liver", "Ercc1(-/D)1 Small intestine", 
                 "Ercc1(-/D)2 Liver", "Ercc1(-/D)2 Small intestine", 
                 "Ercc1(-/D)3 Liver", "Ercc1(-/D)3 Small intestine",  
                 "WT1 Liver", 
                 "WT2 Liver","WT2 Small intestine",
                 "WT3 Liver", "WT3 Small intestine")
vcfs = read_vcfs_as_granges(c(vcf_files,vcf_files_lowVAF), c(sample_names,sample_names_lowVAF), genome = ref_genome)

#2 Check if chromosome names are uniform
all(seqlevels(vcfs[[1]]) %in% seqlevels(get(ref_genome)))

#3 only select autosomal chromosomes
auto = extractSeqlevelsByGroup(species="Mus_musculus", style="UCSC", group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))




# ---- MUTATION ANALYSES: 7 types ----

#1 Get 6 types
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

#2 Plot relative spectrum per 6 mutation types 
plot_spectrum(type_occurrences, by = mousetype, legend = T, CT = T)
p6 <- plot_spectrum(type_occurrences, by = plottype, legend = T, CT = T)
p7 <- plot_spectrum(type_occurrences, by = type, legend = T, CT = T)
#ggsave(paste(outdir,"spectrum_6type_lowVAF.pdf",sep=""), p6, width = 12, height = 10)
#ggsave(paste(outdir,"spectrum_6type_lowVAF_pertype.pdf",sep=""), p7, width = 12, height = 10)

#4 COllapse muts
collapsed_spectrum <- data.frame(colSums(type_occurrences[c(7,8,10),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(1,3,5),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(9,11),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(2,4,6),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(18,19,21),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(12,14,16),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(20,22),], na.rm = FALSE, dims = 1),
                                 colSums(type_occurrences[c(13,15,17),], na.rm = FALSE, dims = 1))
colnames(collapsed_spectrum) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine",
                                  "WT Liver lowVAF","Ercc1(-/D) Liver lowVAF","WT Small intestine lowVAF","Ercc1(-/D) Small intestine lowVAF")

#6 Cosine similarity
cosine_similarity_spectrum <- as.matrix(cos_sim_matrix(collapsed_spectrum[-3,],collapsed_spectrum[-3,]))
p8 <- plot_cosine_heatmap(cosine_similarity_spectrum,cluster_rows = F, plot_values = T)
#ggsave(paste(outdir,"cossim_spectrum_6type_lowVAF.pdf",sep=""), p8, width = 6, height = 5)

#7 Chisquares
#WT liver
p_chi_1 <- chisq.test(collapsed_spectrum[,c(1,5)])$p.value
#MUTLIVER
p_chi_2 <- chisq.test(collapsed_spectrum[,c(2,6)])$p.value
#WT si
p_chi_3 <- chisq.test(collapsed_spectrum[,c(3,7)])$p.value
#MUT si
p_chi_4 <- chisq.test(collapsed_spectrum[,c(4,8)])$p.value

#3 Adjust for multiple testing
p.adjust(c(p_chi_1,p_chi_2,p_chi_3,p_chi_4), method = "fdr")




# ---- MUTATIONAL PROFILES: 96 types ----

#1 Get 96 mutation types
mut_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

#2 COllapse muts
collapsed_mutmatrix <- data.frame(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(12:22)], na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine","lowVAF")

#3 plot 96 types
p1 <- plot_96_profile(collapsed_mutmatrix, ymax = 0.2, condensed = TRUE)
#ggsave(paste(outdir,"spectrum_96type_lowVAF.pdf",sep=""), p1, width = 12, height = 10)

#4 Cosine similarity
cosine_similarity <- as.matrix(cos_sim_matrix(collapsed_mutmatrix,collapsed_mutmatrix))
p5 <- plot_cosine_heatmap(cosine_similarity,cluster_rows = F, plot_values = T)
#ggsave(paste(outdir,"cossim_spectrum_96type_lowVAF.pdf",sep=""), p5, width = 5, height = 4)

#5 Difference in 96 types
plot_compare_profiles(data.frame(rowSums(collapsed_mutmatrix[,c(1:4)], na.rm = FALSE, dims = 1))[,1],
                      collapsed_mutmatrix[,5],
                      profile_names = c("SNV","low VAF"))

#6 Contribution cosmic signatures
cancer_signatures = read.table(paste(indir,"SNV/cancersignatures.txt",sep=""), sep="\t")
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
cancer_signatures = as.matrix(cancer_signatures[,4:33])
fit_res_cancersigs = fit_to_signatures(collapsed_mutmatrix, cancer_signatures)
impt.sigs <- c(1,3,8,9,10,11,14,18,29,30)
select = which(rowSums(fit_res_cancersigs$contribution) > 100)
# > 100 = same signatures as the ones used

#7 plot contribution
p2 <- plot_contribution(fit_res_cancersigs$contribution[impt.sigs,], coord_flip = T, signatures = cancer_signatures[,impt.sigs], mode = "absolute")
p3 <- plot_contribution(fit_res_cancersigs$contribution[impt.sigs,], coord_flip = T, signatures = cancer_signatures[,impt.sigs], mode = "relative")
#ggsave(paste(outdir,"refit_importantcosmic_lowVAF.pdf",sep=""), p2, width = 8, height = 5)
#ggsave(paste(outdir,"refit_importantcosmic_lowVAF_relative.pdf",sep=""), p3, width = 8, height = 5)

#8 cossim to signatures
cossim_cancer_signatures <- cos_sim_matrix(collapsed_mutmatrix, cancer_signatures)
clustered_cancer_signatures <- cluster_signatures(cancer_signatures)
sig_order = colnames(cancer_signatures)[clustered_cancer_signatures$order]
p4 <- plot_cosine_heatmap(cossim_cancer_signatures,cluster_rows = T,plot_values = T,col_order = sig_order)
#ggsave(paste(outdir,"cosine_allcosmic_lowVAF.pdf",sep=""), p4, width = 20, height = 5)

#9 Signature contribution per sample type
collapsed_mutmatrix_pertype <- data.frame(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(18,19,21)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(12,14,16)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(20,22)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(13,15,17)], na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix_pertype) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine",
                                           "WT Liver lowVAF","Ercc1(-/D) Liver lowVAF","WT Small intestine lowVAF","Ercc1(-/D) Small intestine lowVAF")
fit_res_cancersigs_pertype = fit_to_signatures(collapsed_mutmatrix_pertype, cancer_signatures)
p9 <- plot_contribution(fit_res_cancersigs_pertype$contribution[impt.sigs,], coord_flip = T, signatures = cancer_signatures[,impt.sigs], mode = "absolute")
#ggsave(paste(outdir,"refit_importantcosmic_lowVAF_pertype.pdf",sep=""), p9, width = 8, height = 5)
