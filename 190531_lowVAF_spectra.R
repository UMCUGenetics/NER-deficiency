# @Date: 27 march 2019
# @Author: Myrthe Jager
# @Modified: 31 may 2019
# @Description: Plot mutational spectra of variants with VAF < 0.3
# Abbreviations: 




# ---- GET STARTED ----

#1 Load required packages
library(MutationalPatterns)
library("ggplot2")

#2 Define input and output directory
dir = ""
indir = paste(dir,"ERCC1/Data/",sep="")
outdir = paste(dir,"ERCC1/Results/SNV/VAFbelow03/",sep="")

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
p1.1 <- plot_spectrum(type_occurrences, by = type, legend = T, CT = T)
#ggsave(paste(outdir,"spectrum_6type_lowVAF_pertype.pdf",sep=""), p1.1, width = 12, height = 10)





# ---- MUTATIONAL PROFILES: 96 types ----

#1 Get 96 mutation types
mut_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

#2 Collapse muts
collapsed_mutmatrix <- data.frame(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(18,19,21)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(12,14,16)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(20,22)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(13,15,17)], na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine",
                                   "WT Liver lowVAF","Ercc1(-/D) Liver lowVAF","WT Small intestine lowVAF","Ercc1(-/D) Small intestine lowVAF")

#3 plot 96 types
p1.2 <- plot_96_profile(collapsed_mutmatrix, ymax = 0.1, condensed = TRUE)
#ggsave(paste(outdir,"spectrum_96type_lowVAF.pdf",sep=""), p1.2, width = 12, height = 18)