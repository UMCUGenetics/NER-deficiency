# @Date: 18 January 2018
# @Modified: 12 March 2018
# @Author: Myrthe Jager
# @Description: Analysis mutation spectra of old mice (small intestine): contribution signature 8 and similarity to NMF (from ERCC1) sigs
# Abbreviations: SI = Small intestine; SNV = single nucleotide variant (point mutation), WT = wildtype, MUT = mutant (delta min)




# ---- GET STARTED ----

#1 Install & load required packages
#biocLite("MutationalPatterns")
library(MutationalPatterns)
#install.packages("ggplot2")
library("ggplot2")
library("gridExtra")
library(grid)

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/Oldmice/"

#3 Functions
#-

#4 Install and load mouse reference genome
ref_genome_old = "BSgenome.Mmusculus.UCSC.mm9"
#biocLite(ref_genome_old)
library(ref_genome_old, character.only = T)
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
#biocLite(ref_genome)
library(ref_genome, character.only = T)




# ---- GET VCFS OF OLD MICE: SI data ----

#1 Get VCFS
vcf_files_old = list.files(paste(indir,"Oldmice/",sep=""), full.names = T)
sample_names_old = c("MOUSEA568-A5", "MOUSEA568-A6", 
                 "MOUSEA568-A8", "MOUSELGR123-LGR1", 
                 "MOUSELGR123-LGR2", "MOUSELGR123-LGR3")
vcfs_old = read_vcfs_as_granges(vcf_files_old, sample_names_old, genome = ref_genome_old)

#2 Check if chromosome names are uniform
all(seqlevels(vcfs_old[[1]]) %in% seqlevels(get(ref_genome_old)))

#3 only select autosomal chromosomes
auto = extractSeqlevelsByGroup(species="Mus_musculus", style="UCSC", group="auto")
vcfs_old = lapply(vcfs_old, function(x) keepSeqlevels(x, auto))




# ------ GET VCFS of ERCC1 WT and MUT ------

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
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))




# ---- MUTATION TYPES ----

#1 Get 6 types
type_occurrences_old = mut_type_occurrences(vcfs_old, ref_genome_old)
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

#2 Get 96 types
mut_matrix_old = mut_matrix(vcf_list = vcfs_old, ref_genome = ref_genome_old)
mut_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

#3 Collapse
collapsed_mutmatrix_old <- data.frame(rowSums(mut_matrix_old, na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix_old) <- "oldmice_smallintestine"




# ---- SIG 8 CONTRIBUTION IN WT SI ----

#1 Retrieve signatures from pan-cancer study Alexandrov et al.
cancer_signatures = read.table(paste(indir,"SNV/cancersignatures.txt",sep=""), sep="\t")

#2 Reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]

#3 Only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

#4 Refit
fit_res_old = fit_to_signatures(mut_matrix_old, cancer_signatures[,c(1,8,10,14,18,29,30)])
fit_res_new = fit_to_signatures(mut_matrix, cancer_signatures[,c(1,8,10,14,18,29,30)])

#5 Determine Sig 8 contribution
sig8_old <- data.frame(t(t(fit_res_old$contribution[2,])),c(98,98,98,116,116,116))
sig8_new <- data.frame(t(t(fit_res_new$contribution[2,])),rep(16,11))
colnames(sig8_old) <- c("mutations","age")
colnames(sig8_new) <- c("mutations","age")
sig8_old$perweek <- as.numeric(sig8_old[,1])/as.numeric(sig8_old[,2])
sig8_new$perweek <- as.numeric(sig8_new[,1])/as.numeric(sig8_new[,2])
sig8_rate_df<- rbind(sig8_old,sig8_new[c(9,11),])

#5 Plot sig 8 contribution
sig8plot_old <- ggplot(sig8_rate_df, aes(x = age,y=mutations)) +
  geom_point(shape=1) +
#  geom_smooth(method=lm) + 
  xlim(0,150) +
  ylim(0,400)+
  xlab("Age (weeks)") +
  ylab("Signature 8 mutations") +
  ggtitle("Number of Signature 8 mutations\naccumulated in ASCs\nof mouse small intestine") +
  theme(plot.title = element_text(size = rel(1.5),lineheight=.8,hjust = 0.5))
#ggsave(paste(outdir,"sig8_absolutecontribution_witholdmice.pdf",sep=""),plot = sig8plot_old, width = 4, height = 4)




# ---- COSINE SIMILARITY: NMF SIGS ERCC1 ----

#1 Retrieve Signatures NMF ERCC1
nmf_sigs <- read.table(paste(indir,"SNV/NMF_2signatures.txt",sep=""), sep = "\t")
#nmf_sigs <- nmf_res2$signatures

#2 Generate one big df
df.big <- cbind(nmf_sigs,collapsed_mutmatrix_old,cancer_signatures)

#3 Calculate cosine similarity
cosine_similarity_old <- cos_sim_matrix(df.big,df.big)

#4 Plot cosine similarity
cosine_plot_old <- plot_cosine_heatmap(cosine_similarity_old,cluster_rows = T, plot_values = T)

#5 Save plot
#ggsave(paste(outdir,"NMF_2signatures_cossim_witholdmice.pdf",sep=""), plot = cosine_plot_old, height = 16, width = 16)
