# @Date: 25 July 2017
# @Author: Myrthe Jager
# @Modified: 13 March 2018, 27 March 2018
# @Description: Genome-wide mutation pattern analysis ERCC1 delta min mouse organoids
# @Type: MUT & WT
# @Tissue: liver & small intestine
# Abbreviations: MUT/D- = Ercc1 delta min; WT = Wild-type, SI = Small intestine; SNV = single nucleotide variant (point mutation)




# ---- GET STARTED ----

#1 Install & load required packages
#install.packages("devtools")
library(devtools)
#options(unzip = 'internal')
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("BiocInstaller")
library(BiocInstaller)
#biocLite("MutationalPatterns")
#biocLite("GenomeInfoDbData")
#biocLite("GenomeInfoDb")
#install_github("UMCUgenetics/MutationalPatterns",ref = "develop")
library(MutationalPatterns)
library("gridExtra")
library(grid)
#install.packages("https://cran.r-project.org/src/contrib/ggplot2_2.2.1.tar.gz", repo=NULL, type="source")
library("ggplot2")
#install.packages("https://cran.r-project.org/web/packages/R.matlab/index.html")
library(grDevices)
library(reshape2)
library(NMF)
#biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
#biocLite("biomaRt")
library(biomaRt)

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/SNV/"

#3 Functions
## Function 1: PER TYPE: STATISTICS ##
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Function 2: PLOT ABSOLUTE MUTATIONAL SPECTRA (96 TYPES) ##
## Adjusted from Mutational Patterns package
plot_96_profile_absolute = function(mut_matrix, colors, ymax = 0.2, condensed = FALSE) {
  # Check color vector length
  # Colors for plotting
  COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")
  
  COLORS7 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#E98C7B", "#D4D2D2", "#ADCC54",
    "#F0D0CE")
  
  SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
  SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
  SUBSTITUTIONS_192 = rep(SUBSTITUTIONS, each=32)
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  
  # combine substitutions and context in one 
  TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")
  
  STRAND = rep(c("U","T"), 96)
  DNA_BASES = c("A", "C", "G", "T")
  
  if(missing(colors)){colors=COLORS6}
  if(length(colors) != 6){stop("Provide colors vector with length 6")}
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each=16)
  
  # Replace mutated base with dot to get context
  substring(context, 2, 2) = "."
  
  # Construct dataframe
  df = data.frame(substitution = substitution, context = context)
  rownames(mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  
  plot = ggplot(data=df3, aes(x=context,
                                y=value,
                                fill=substitution,
                                width=1)) +
      geom_bar(stat="identity", colour="black", size=.2) +
      scale_fill_manual(values=colors) +
      facet_grid(variable ~ substitution) +
      ylab("Absolute contribution") +
      coord_cartesian(ylim=c(0,ymax)) +
      scale_y_continuous(breaks=seq(0, ymax, 10)) +
      # no legend
      guides(fill=FALSE) +
      # white background
      theme_bw() +
      # format text
      theme(axis.title.y=element_text(size=12,vjust=1),
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=12),
            axis.text.x=element_text(size=5,angle=90,vjust=0.4),
            strip.text.x=element_text(size=9),
            strip.text.y=element_text(size=9),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.spacing.x = unit(0, "lines"))
  return(plot)
}

## Function 3: Rename CHROM to UCSC style ##
rename_chrom = function(granges, style = "UCSC")
{
  # rename mitochondrial DNA manually
  seqlevels(granges)[seqlevels(granges)=="chrMT"] = "chrM"
  
  # get chromosome style
  chrom_style = mapSeqlevels(seqlevels(granges), style)
  
  # removing NA cases
  chrom_style = chrom_style[complete.cases(chrom_style)] 
  
  # rename chromosome names (seqlevels)
  res = renameSeqlevels(granges, chrom_style)
  
  return(res)
}

## Function 4: Convert bed to Granges ##
bed_to_granges = function(bed_files, region_names) 
{ 
  if (length(bed_files) != length(region_names)) 
    stop("Provide the same number of names as bed files") 
  
  granges_list = list() 
  for(i in 1:length(bed_files)) 
  { 
    bed_file = bed_files[i] 
    bed = read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
    chr = paste("chr", bed[,1], sep="") 
    
    # Convert BED (0-based) start postion to Granges (1-based) 
    start = bed[,2] + 1 
    
    # In BED end position is excluded, in Granges end position is 
    # included -> +1 -1 -> no conversion needed 
    end = bed[,3] 
    
    new_bed = GRanges(chr, IRanges(start,end)) 
    new_bed = list(new_bed) 
    names(new_bed) = region_names[i] 
    granges_list = c(granges_list, new_bed) 
  } 
  
  return(granges_list) 
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




# ------ MUTATION NUMBER ------

#1 Load corrected number file
snvcorrectedmutnumber <- read.delim2(paste(indir,"callable_genome/Callable_ERCC1.txt",sep=""))
snvcorrectedmutnumber$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(snvcorrectedmutnumber$Surveyed.percentage)))
snvcorrectedmutnumber$type <- c(rep(c("Ercc1-/D Liver","Ercc1-/D Small intestine"),3),rep(c("WT Liver","WT Small intestine"),3))

#2 Add double nucleotide counts (later on: subtract these from SNVs)
#2A Import table: IGV-verfied double base substitutions
checkeddouble <- read.delim2(paste(indir,"SNV/dinucleotides_igv.txt",sep=""))
checkeddouble <- checkeddouble[which(checkeddouble$eyeball == "TP"),]
#2B Count double base substitutions & add to callable genome table
snvcorrectedmutnumber$uncorr.double <- c(as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)1 Liver"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)1 Small intestine"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)2 Liver"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)2 Small intestine"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)3 Liver"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)3 Small intestine"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "WT1 Liver"))),
                                            NA,
                                            as.numeric(length(which(checkeddouble$mousetype == "WT2 Liver"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "WT2 Small intestine"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "WT3 Liver"))),
                                            as.numeric(length(which(checkeddouble$mousetype == "WT3 Small intestine"))))
remove(checkeddouble)

#3 Extract SNV numbers from file
snvcorrectedmutnumber$uncorr.snv.di <- NA
for (i in vcf_files) {
  df <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
  snv_number <- nrow(df)
  
  #3B Add to the surveyed file
  snvcorrectedmutnumber[which(as.factor(toupper(substr(snvcorrectedmutnumber$Sample,1,5))) == toupper(substr(sample,6,10))),]$uncorr.snv.di <- snv_number
}
remove(i,sample,df,snv_number)

#4 Subtract Dinucleotide changes (these are in the SNV vcfs; but should only be counted as double nucleotide substitution, not both)
snvcorrectedmutnumber$uncorr.snv <- snvcorrectedmutnumber$uncorr.snv.di-(snvcorrectedmutnumber$uncorr.double *2)

#5 Correct mutationnumber for surveyed area (extrapolate SNVs to autosomal genome)
snvcorrectedmutnumber$corr.snv <- (snvcorrectedmutnumber$uncorr.snv/snvcorrectedmutnumber$Surveyed.percentage)*100

#6 Calculate mutation rate per week
snvcorrectedmutnumber$snvrate <- snvcorrectedmutnumber$corr.snv/16
#write.table(snvcorrectedmutnumber,file = paste(outdir,"snvnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)

#7 Generate data frame with all 12 samples (Ercc1/WT & Liver/SI)
mutnumber <- data.frame(
  type = factor(c("WT Liver","Ercc1-/D Liver","WT Small intestine","Ercc1-/D Small intestine"),
                levels=c("WT Liver","Ercc1-/D Liver","WT Small intestine", "Ercc1-/D Small intestine")),
  mouse = c("WT1 Liver","Ercc1(-/D)1 Liver","WT1 Small intestine","Ercc1(-/D)1 Small intestine",
            "WT2 Liver","Ercc1(-/D)2 Liver","WT2 Small intestine","Ercc1(-/D)2 Small intestine",
            "WT3 Liver","Ercc1(-/D)3 Liver","WT3 Small intestine","Ercc1(-/D)3 Small intestine"), 
  mutations = snvcorrectedmutnumber[c(7,1,8,2,9,3,10,4,11,5,12,6),]$snvrate
)




# ---- STATISTICAL ANALYSIS: SNV rate ----

#1 Welch Two Sample t-test
pliver.snv<-t.test(mutnumber[grep("Liver",mutnumber$type),"mutations"] ~ mutnumber[grep("Liver",mutnumber$type),"type"])$p.value # p-value = 0.001164
psi.snv<-t.test(mutnumber[grep("Small",mutnumber$type),"mutations"] ~ mutnumber[grep("Small",mutnumber$type),"type"])$p.value # NS

#2 Adjust for multiple testing
p.adjust(c(pliver.snv,psi.snv), method = "fdr") #0.002644993 0.262298805





# ---- PLOT: SNV rate ----

#1 Per ASC type: get statistics (for plotting)
muts <- summarySE(mutnumber[-3,], measurevar="mutations",groupvars="type")

#2 Per ASC type (4 types: Ercc1/WT & Liver/SI)
#2A Plot
mutnumberplot <- ggplot(data=muts, aes(x=type, y=mutations, fill= type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin=mutations-sd, ymax=mutations+sd),
                width=.07,                    # Width of the error bars
                position=position_dodge(.9)) +
  xlab(" Ercc1") + 
  ylab("Base substitutions\nper autosomal genome per week\n") +
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=50) +
  scale_fill_manual(values=c("lightcoral","skyblue3",
                             "lightcoral","skyblue3")) +
  guides(fill=FALSE) +  #geen legenda
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())
#2B Save plot
#ggsave(paste(outdir,"SNVnumber.pdf",sep=""), mutnumberplot, width = 8, height = 8)

#3 Per sample (11 samples)
#3A Plot
mutnumberpersampleplot <- ggplot(data=mutnumber, aes(x=type, y=mutations, fill= mouse)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab(" Ercc1") + 
  ylab("Base substitutions per autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("skyblue3","skyblue3","lightskyblue2","lightskyblue2","steelblue3","steelblue3",
                             "lightpink","lightpink","lightcoral","lightcoral","indianred3","indianred3")) +
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  expand_limits(y=50) +
  guides(fill=FALSE) +  #geen legenda
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())
#3B Save
#ggsave(paste(outdir,"SNVnumber_permouse.pdf",sep=""), mutnumberpersampleplot, width = 8, height = 8)




# ---- MUTATION SPECTRA ----

#1 Get 6 types
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

#2 Plot relative spectrum per 6 mutation types 
#2A Plot (per mousetype with C>T at CpG separate and per animal)
spectrum6typeplotwithCT = plot_spectrum(type_occurrences, by = mousetype, legend = T, CT = T)
spectrum6typeplotpermouse = plot_spectrum(type_occurrences, by = sample_names, legend = T, CT = T)
#2B Save plots
#ggsave(paste(outdir,"spectrum_6type_withCTatCpG.pdf",sep=""), spectrum6typeplotwithCT, width = 8, height = 8)
#ggsave(paste(outdir,"extra/spectrum_6type_permouse.pdf",sep=""), spectrum6typeplotpermouse, width = 10, height = 10)

#3 Get 96 mutation types
mut_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

#4 Absolute spectra
#4A Plot
liverd96 = plot_96_profile_absolute(mut_matrix[,c(1,3,5)], ymax = 50)
liverwt96 = plot_96_profile_absolute(mut_matrix[,c(7,8,10)], ymax = 50)
sid96= plot_96_profile_absolute(mut_matrix[,c(2,4,6)], ymax = 50)
siwt96 = plot_96_profile_absolute(mut_matrix[,c(9,11)], ymax = 50)
#4B Save plot
#ggsave(paste(outdir,"spectrum_96type_permouse_abs.pdf",sep=""), grid.arrange(liverd96, sid96, liverwt96,siwt96, nrow=2, ncol =2), width = 12, height = 12)

#5 Relative spectra
#5A Plot
liverd96_rel = plot_96_profile(mut_matrix[,c(1,3,5)], ymax = 0.1,condensed = T)
liverwt96_rel = plot_96_profile(mut_matrix[,c(7,8,10)], ymax = 0.1,condensed = T)
sid96_rel= plot_96_profile(mut_matrix[,c(2,4,6)], ymax = 0.1,condensed = T)
siwt96_rel = plot_96_profile(mut_matrix[,c(9,11)], ymax = 0.1,condensed = T)
#5B Save plot
#ggsave(paste(outdir,"extra/spectrum_96type_permouse_rel.pdf",sep=""), grid.arrange(liverd96_rel, sid96_rel, liverwt96_rel,siwt96_rel, nrow=2, ncol =2), width = 12, height = 12)




# ---- STATISTICAL ANALYSIS: Mutation spectra----

#1 Generate data frames for statistical testing: 7 types
#new_typeoccurences <- t(type_occurrences[,1:6])
new_typeoccurences <- t(type_occurrences[,c(1:2,4:8)])
liver.spectrum <- data.frame(wtliv = rowSums(new_typeoccurences[,c(7,8,10)]),
                             mutliv = rowSums(new_typeoccurences[,c(1,3,5)]))
si.spectrum <- data.frame(wtsi = rowSums(new_typeoccurences[,c(9,11)]),
                          mutsi = rowSums(new_typeoccurences[,c(2,4,6)]))

#2 chi square 6 types
p_chi_liver_6type <- chisq.test(liver.spectrum)$p.value
p_chi_si_6type <- chisq.test(si.spectrum)$p.value

#3 Adjust for multiple testing
p.adjust(c(p_chi_liver_6type,p_chi_si_6type), method = "fdr") #0.040948305 0.004717594




# ---- COSINE SIMILARITY OF Mutational spectra ----

#1 Generate df for calculating 
#1A Collapse mutmatrix
collapsed_mutmatrix <- data.frame(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1),
                                  rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1))
colnames(collapsed_mutmatrix) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")
#1B Correct for number of samples (n=2 WT SI, others n=3)
collapsed_mutmatrix_combined <- data.frame(round(rowSums(mut_matrix[,c(7,8,10)], na.rm = FALSE, dims = 1)/3),
                                           round(rowSums(mut_matrix[,c(1,3,5)], na.rm = FALSE, dims = 1)/3),
                                           round(rowSums(mut_matrix[,c(9,11)], na.rm = FALSE, dims = 1)/2),
                                           round(rowSums(mut_matrix[,c(2,4,6)], na.rm = FALSE, dims = 1)/3))
colnames(collapsed_mutmatrix_combined) <- c("WT Liver","Ercc1(-/D) Liver","WT Small intestine","Ercc1(-/D) Small intestine")

#2 plot 96 profile per type
#ggsave(paste(outdir,"extra/spectrum_96type_pertype_rel.pdf",sep=""), plot = plot_96_profile(collapsed_mutmatrix,condensed = T), width = 8, height = 8)
#ggsave(paste(outdir,"extra/spectrum_96type_pertype_abs.pdf",sep=""), plot = plot_96_profile_absolute(collapsed_mutmatrix_combined,ymax=50), width = 8, height = 8)

#3 Calculate cosine similarity
#3A Combined spectra per ASCtype (4 types: WT/D- and Liver/SI)
cosine_similarity <- as.matrix(cos_sim_matrix(collapsed_mutmatrix,collapsed_mutmatrix))
#3B Per mouse
mut_matrix = as.matrix(mut_matrix) 
mut_matrix_ordered = mut_matrix[,c(7,8,10,1,3,5,9,11,2,4,6)]
cosine_similarity_permouse <- as.matrix(cos_sim_matrix(mut_matrix_ordered,mut_matrix_ordered))
#write.table(mut_matrix_ordered,file =paste(outdir,"mutational_profile.txt",sep=""),sep ="\t", row.names = F, col.names = T)

#4 plot cosine similarity
#4A Per ASCtype (4 types: WT/D- and Liver/SI)
cossimplot <- plot_cosine_heatmap(cosine_similarity,cluster_rows = F, plot_values = T)
#4B Per mouse
cossimplot_permouse <- plot_cosine_heatmap(cosine_similarity_permouse,cluster_rows = F, plot_values = F)
#4C Save plots
#ggsave(paste(outdir,"cosinesimilarity_mutationspectra.pdf",sep=""), plot = cossimplot, width = 6, height = 6)
#ggsave(paste(outdir,"cosinesimilarity_mutationspectra_permouse.pdf",sep=""), plot = cossimplot_permouse, width = 6, height = 6)




# ---- COSMIC SIGNATURES: COSINE SIMILARITY ----

#1 Generate table with cosmic signatures
#1A Retrieve signatures from pan-cancer study Alexandrov et al.
#url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
#cancer_signatures = read.table(url, sep = "\t", header = T)
#write.table(cancer_signatures, file  = paste(indir,"SNV/cancersignatures.txt",sep=""), sep = "\t") # v. 13 March 2018
cancer_signatures = read.table(paste(indir,"SNV/cancersignatures.txt",sep=""), sep="\t")
#1B Reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
#1C Only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

#2 Calculate cosine similarity to COSMIC cancer signatures
#2A For collapsed per mousetype
cossim_cancer_signatures <- cos_sim_matrix(collapsed_mutmatrix, cancer_signatures)
#2B For each ASC separately
cossim_cancer_signatures_permouse <- cos_sim_matrix(mut_matrix_ordered, cancer_signatures)

#3 Plot
#3A Order sigs for plot
clustered_cancer_signatures <- cluster_signatures(cancer_signatures)
sig_order = colnames(cancer_signatures)[clustered_cancer_signatures$order]
#3B Plot cossims
cossimcosmicplot <- plot_cosine_heatmap(cossim_cancer_signatures,cluster_rows = T,plot_values = T,col_order = sig_order)
cosine_plot <- plot_cosine_heatmap(cossim_cancer_signatures_permouse, cluster_rows = F, col_order = sig_order,plot_values = F)
#3C Save plots
#ggsave(paste(outdir,"cosinesimilarity_tocosmicsignatures.pdf",sep=""), plot = cossimcosmicplot, width = 15, height = 4)
#ggsave(paste(outdir,"cosinesimilarity_tocosmicsignatures_permouse.pdf",sep=""), plot = cosine_plot, width = 15, height = 8)




# ---- COSMIC SIGNATURES: CONTRIBUTION ----

#1 Refit
fit_res_cancersigs = fit_to_signatures(mut_matrix_ordered, cancer_signatures)
fit_res_cancersigs.collapsed = fit_to_signatures(collapsed_mutmatrix, cancer_signatures)
fit_res_collapsed_combined= fit_to_signatures(collapsed_mutmatrix_combined, cancer_signatures)

#3 Calculate relative contribution of signature to mutation spectrum 
#3A of 4 ASC types 
collapsed_mutmatrix.fit_rel <- fit_res_cancersigs.collapsed$contribution
for (i in 1:ncol(collapsed_mutmatrix.fit_rel)){
  collapsed_mutmatrix.fit_rel[,i] <- collapsed_mutmatrix.fit_rel[,i]/sum(collapsed_mutmatrix.fit_rel[,i])
}
remove(i)
#3B per ASC
fit_res_rel_cont <- fit_res_cancersigs$contribution
for (i in 1:ncol(fit_res_rel_cont)){
  fit_res_rel_cont[,i] <- fit_res_rel_cont[,i]/sum(fit_res_rel_cont[,i])
}
remove(i)

#4 Plot relative contribution of sigs
#4A Prepare tables for plot
collapsed_mutmatrix.fit_rel <- as.matrix(t(collapsed_mutmatrix.fit_rel))
#4B Plot
contributioncosmicplot.all <- plot_cosine_heatmap(collapsed_mutmatrix.fit_rel,cluster_rows = T, plot_values = T,col_order = sig_order)
contributioncosmicplot.all.bar <- plot_contribution(t(collapsed_mutmatrix.fit_rel))
contributioncosmicplot.all_notcombined <- plot_cosine_heatmap(t(fit_res_rel_cont),cluster_rows = F, plot_values = F,col_order = sig_order)
#4C Save plots
#ggsave(paste(outdir,"refit_allcosmic_permouse.pdf",sep=""), plot = contributioncosmicplot.all_notcombined, width = 15, height = 8)
#ggsave(paste(outdir,"refit_allcosmic.pdf",sep=""), plot = contributioncosmicplot.all, width = 15, height = 4)
#ggsave(paste(outdir,"refit_allcosmic_relative.pdf",sep=""), plot = contributioncosmicplot.all.bar, width = 8, height = 4)

#5 Cosine sim reconstructed
#5A All mice separate
reconstructed_allcosmicsigs <- fit_res_cancersigs$reconstructed
colnames(reconstructed_allcosmicsigs) <- c("rec WT1 Liver","rec WT2 Liver","rec WT3 Liver",
                                           "rec Ercc1(-/D)1 Liver","rec Ercc1(-/D)2 Liver","rec Ercc1(-/D)3 Liver",
                                           "rec WT2 Small intestine","rec WT3 Small intestine",
                                           "rec Ercc1(-/D)1 Small intestine","rec Ercc1(-/D)2 Small intestine","rec Ercc1(-/D)3 Small intestine")
cossimcosmicplot.all <- plot_cosine_heatmap(cos_sim_matrix(reconstructed_allcosmicsigs,mut_matrix_ordered),plot_values = TRUE, cluster_rows = FALSE)
#ggsave(paste(outdir,"extra/similarity_reconstructed_allcosmic.pdf",sep=""), plot = cossimcosmicplot.all, width = 8, height = 8)
#5B collapsed per type
reconstructed_allcosmicsigs_collapsed <- fit_res_cancersigs.collapsed$reconstructed
colnames(reconstructed_allcosmicsigs_collapsed) <- c("rec WT Liver","rec Ercc1(-/D) Liver","rec WT Small intestine","rec Ercc1(-/D) Small intestine")
cossimcosmicplot.all.collapsed <- plot_cosine_heatmap(cos_sim_matrix(reconstructed_allcosmicsigs_collapsed,collapsed_mutmatrix),plot_values = TRUE, cluster_rows = FALSE)
#ggsave(paste(outdir,"extra/similarity_reconstructed_allcosmic_collapsed.pdf",sep=""), plot = cossimcosmicplot.all.collapsed, width = 4, height = 4)




# ---- Define number of sigs for refitting ----

#1 Generate df
test.sigs <- data.frame(fit_res_collapsed_combined$contribution)
test.sigs$sum <- rowSums(test.sigs)
test.sigs$rank <- NA
test.sigs$sig <- NA

#2 Rank signatures, based on contribution
rank.sigs <- order(-test.sigs$sum)
for(i in rank.sigs) {
  sig <- rank.sigs[i]
  test.sigs[sig,]$sig <- as.numeric(sig)
  
  rank <- i
  test.sigs[sig,]$rank <- rank
  remove(sig,rank)
}
test.sigs <- test.sigs[order(-test.sigs$sum),]

#3 Generate new empty df
df.sigs.test <- data.frame(n.sigs = NA,
                           wt.liver = NA,
                           mut.liver = NA,
                           wt.si = NA,
                           mut.si = NA,
                           mean = NA)

#4 Calculate how much of the original profile is reconstructed by the refitted data
# First 2 sigs, then 3, etc.
for(i in 2:nrow(test.sigs)) {
  sigs.to.test <- test.sigs[1:i,]$sig
  cancer.sigs <- cancer_signatures[,sigs.to.test]
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

#5 Plot
df.sigs.test.forplot <- melt(df.sigs.test[,-6], id.vars = "n.sigs")
plot.sigs <- ggplot(df.sigs.test.forplot, aes(x=n.sigs,y=value,fill= variable)) +
  geom_bar(stat = "identity", position = position_dodge())+ 
  ggtitle("Collapsed & adjusted for n mouse") +
  geom_hline(yintercept = 0.90) +
  scale_x_continuous("Signatures", labels = as.character(2:30), breaks = c(2:30)) +
  coord_cartesian(ylim=c(0.8, 1)) +
  theme_bw()
#ggsave(paste(outdir,"nsigs_collapsed_adjustedforn.pdf",sep=""),plot = plot.sigs, width = 12, height = 4)
# 10 signatures! After that, the reconstructed profiles do not improve further



# ---- IMPORTANT COSMIC SIGNATURES ----

#1 Select signatures which are important
impt.sigs <- sort(test.sigs[1:10,]$sig)

#2 Conclusion: 1, 3, 8, 9, 10, 11, 14, 18, 29, and 30

#3 Relative contribution of important signatures
#3A Refit to important signatures
#3A1 per mouse
fit_res_cancersigs_impt <- fit_to_signatures(mut_matrix_ordered,cancer_signatures[,impt.sigs])
#3A2 Per type
fit_res_cancersigs_impt_sum <- fit_to_signatures(collapsed_mutmatrix,cancer_signatures[,impt.sigs])
#3A3 Per type, corrected for number of mice in each group
fit_res_collapsed_combined_impt <- fit_to_signatures(collapsed_mutmatrix_combined,cancer_signatures[,impt.sigs])
#3B Plot refit (per type, corrected for number of mice in each group)
contributioncosmicplot.important <- plot_cosine_heatmap(collapsed_mutmatrix.fit_rel[,impt.sigs],cluster_rows = F, plot_values = T)
#3C Save plot
#ggsave(paste(outdir,"extra/refit_important.pdf",sep=""), plot = contributioncosmicplot.important, width = 6, height = 4)

#4 Plot absolute contribution of sigs
abscontribution_collapsed_combined <- plot_contribution(fit_res_collapsed_combined_impt$contribution, coord_flip = F, signatures = cancer_signatures[,impt.sigs], mode = "absolute")
relcontribution_collapsed_combined <- plot_contribution(fit_res_collapsed_combined_impt$contribution, coord_flip = F, signatures = cancer_signatures[,impt.sigs], mode = "relative")
#ggsave(paste(outdir,"refit_importantcosmic.pdf",sep=""), plot = abscontribution_collapsed_combined,width = 4, height = 4)
#ggsave(paste(outdir,"refit_importantcosmic_relative.pdf",sep=""), plot = relcontribution_collapsed_combined,width = 4, height = 4)

#5 Cosine similarity of reconstructed of refitted to original
#5A reconstructed vs prior to reconstructing (per mouse, important sigs)
fit_res_cancersigs_impt_reconstructed <- fit_res_cancersigs_impt$reconstructed
colnames(fit_res_cancersigs_impt_reconstructed) <- c("rec_WT1 Liver","rec_WT2 Liver","rec_WT3 Liver",
                                             "rec_Ercc1(-/D)1 Liver","rec_Ercc1(-/D)2 Liver","rec_Ercc1(-/D)3 Liver",
                                             "rec_WT2 Small intestine","rec_WT3 Small intestine",
                                             "rec_Ercc1(-/D)1 Small intestine","rec_Ercc1(-/D)2 Small intestine","rec_Ercc1(-/D)3 Small intestine")
cossim_permouse_reconstructed <- plot_cosine_heatmap(cos_sim_matrix(fit_res_cancersigs_impt_reconstructed,mut_matrix_ordered), cluster_rows = F, plot_values = T)
#ggsave(paste(outdir,"extra/similarity_reconstructed_permouse_importantcosmicsigs_original_permouse.pdf",sep=""), plot = cossim_permouse_reconstructed, height = 8, width = 8)
#5B reconstructed collapsed vs prior to reconstructing per mouse (important sigs)
t1 <- fit_res_collapsed_combined_impt$reconstructed
rownames(t1) <- row.names(collapsed_mutmatrix_combined)
colnames(t1) <-c("reconstructed WTliver","reconstructed Ercc1(-/D)Liver","reconstructed WT SI","reconstructed Ercc1(-/D)SI")
cossim_reconstructed_permouse <- plot_cosine_heatmap(cos_sim_matrix(t1,mut_matrix_ordered),plot_values = T, cluster_rows = F)
#ggsave(paste(outdir,"extra/similarity_reconstructed_importantcosmicsigs_original_permouse.pdf",sep=""), plot = cossim_reconstructed_permouse, height = 4, width = 8)
#5C collapsed reconstructued vs collapsed prior to reconstructing (collapsed per type, important sigs, corrected for number of mice per type)
cossim_reconstructed <- plot_cosine_heatmap(cos_sim_matrix(t1,collapsed_mutmatrix_combined),plot_values = T, cluster_rows = F)
#ggsave(paste(outdir,"extra/similarity_reconstructed_importantcosmicsigs_original.pdf",sep=""), plot = cossim_reconstructed, height = 8, width = 8)

#6 Cossim to important signatures
#6A Calculate cossim
cossim_cancer_signatures_imp <- cos_sim_matrix(collapsed_mutmatrix,cancer_signatures[,impt.sigs])
#6B Plot cossim
cosine_plot_imp <- plot_cosine_heatmap(cossim_cancer_signatures_imp,cluster_rows = F, plot_values = T)
#6C Save plot
#ggsave(paste(outdir,"cosinesimilarity_tocosmicsignatures_important.pdf",sep=""), plot = cosine_plot_imp, width = 6, height = 4)




# ---- NMF ----

#1 Define number of signatures
mut_matrix_ordered_2 = mut_matrix_ordered + 0.0001 
#estim.r = NMF::nmf(mut_matrix_ordered_2, rank = 2:10, method = "brunet", nrun = 100, seed = 123456)
#ggsave(paste(outdir,"extra/NMF_estimroutput.pdf",sep=""), plot = plot(estim.r), width = 8, height = 8)

#2 Extract two signatures
nmf_res2 = extract_signatures(mut_matrix_ordered_2, rank = 2)
colnames(nmf_res2$signatures) = c("Signature A", "Signature B")
rownames(nmf_res2$contribution) = c("Signature A", "Signature B")

#3 Plot
profile962 = plot_96_profile(nmf_res2$signatures,ymax = 0.1,condensed = T)
contribution962 = plot_contribution(nmf_res2$contribution, signatures = nmf_res2$signatures, coord_flip = T, mode = "absolute")
#3B Save plot
#ggsave(paste(outdir,"NMF_2signatures.pdf",sep=""), plot = grid.arrange(profile962, contribution962, nrow=1, ncol =2), width = 16, height = 8)

#4 Generate and save table with relative signature contributions
#4A Generate table
relative_nmf <- nmf_res2$signatures
relative_nmf <- rbind(relative_nmf,colSums(relative_nmf))
for(i in 1:nrow(relative_nmf)) {
  relative_nmf[i,1] <- relative_nmf[i,1]/relative_nmf[97,1]
  relative_nmf[i,2] <- relative_nmf[i,2]/relative_nmf[97,2]
}
relative_nmf <- relative_nmf[-97,]
#4B Save table
#write.table(relative_nmf, file  = paste(outdir,"NMF_2signatures.txt",sep=""), sep = "\t")
#write.table(relative_nmf, file  = paste(indir,"SNV/NMF_2signatures.txt",sep=""), sep = "\t")

#5 Cosine similarity reconstructed vs original mutation spectra
t1_2sigs <- nmf_res2$reconstructed
colnames(t1_2sigs) <- c("reconstructed WT1 Liver", "reconstructed WT2 Liver","reconstructed WT3 Liver","reconstructed Ercc1(-/D)1 Liver","reconstructed Ercc1(-/D)2 Liver","reconstructed Ercc1(-/D)3 Liver","reconstructed WT2 Small intestine","reconstructed WT3 Small intestine","reconstructed Ercc1(-/D)1 Small intestine", "reconstructed Ercc1(-/D)2 Small intestine", "reconstructed Ercc1(-/D)3 Small intestine")
cossim_2sigs <- plot_cosine_heatmap(cos_sim_matrix(t1_2sigs,mut_matrix_ordered),plot_values = T, cluster_rows = F)
#ggsave(paste(outdir,"extra/similarity_reconstructed2sigs_original.pdf",sep=""), plot = cossim_2sigs, height = 8, width = 8)
remove(nmf_res2)

#6 Similarity ERCC1 sigs to cosmic sigs
#6A Get signatues
ercc1_signatures = read.table(file = paste(indir,"SNV/NMF_2signatures.txt",sep=""), sep = "\t")
ercc1_signatures <- as.matrix(ercc1_signatures)
#6B Calculate similarity to COSMIC signatures
cancer_signatures_forcossim <- cancer_signatures
rownames(cancer_signatures_forcossim) <- rownames(mut_matrix)
cosine_similarity_nmf <- as.matrix(cos_sim_matrix(ercc1_signatures,cancer_signatures_forcossim))
#6C plot cosine similarity
cosplot_2sigs <- plot_cosine_heatmap(cosine_similarity_nmf,cluster_rows = F, plot_values = T, col_order = sig_order)
#ggsave(paste(outdir,"NMF_2signatures_cossimtocosmic.pdf",sep=""), plot = cosplot_2sigs, height = 4, width = 16)

#7 Refit collapsed to signatures
ercc1_fitted_to_ercc1 <- fit_to_signatures(collapsed_mutmatrix_combined, ercc1_signatures)
contribution.plot.ercc1fitted <- plot_contribution(ercc1_fitted_to_ercc1$contribution, coord_flip = F, signatures = ercc1_signatures, mode = "absolute")
#ggsave(paste(outdir,"NMF_2signatures_fitted.pdf",sep=""), plot = contribution.plot.ercc1fitted, height = 6, width = 4)




# ---- Signature 8 ----

#1 Plot Signature 8
cancer_plot<-data.frame(cancer_signatures[,8],ercc1_signatures[,1])
colnames(cancer_plot) <- c("Signature 8","Signature 8*")
#ggsave(paste(outdir,"signature8.pdf",sep=""), plot = plot_96_profile(cancer_plot,ymax = 0.05, condensed = T))

sig8plot <- plot_compare_profiles(profile1 = cancer_plot[,1],profile2 = cancer_plot[,2],profile_names = c("Signature 8","Signature 8*"), profile_ymax = 0.1, condensed = TRUE,diff_ylim = c(-0.05, 0.05))
#ggsave(paste(outdir,"signature8_and8like.pdf",sep=""), plot = sig8plot, height = 6, width = 8)




# ---- RAINFALL PLOT ----

#1 define chromosomes of interest
chromosomes = seqnames(get(ref_genome))[1:19]

#2 Plot rainfall
rainmut1liver = plot_rainfall(vcfs[[1]], title = names(vcfs[1]), chromosomes = chromosomes, cex = 1)
rainmut2liver = plot_rainfall(vcfs[[3]], title = names(vcfs[3]), chromosomes = chromosomes, cex = 1)
rainmut3liver = plot_rainfall(vcfs[[5]], title = names(vcfs[5]), chromosomes = chromosomes, cex = 1)
rainmut1si = plot_rainfall(vcfs[[2]], title = names(vcfs[2]), chromosomes = chromosomes, cex = 1)
rainmut2si = plot_rainfall(vcfs[[4]], title = names(vcfs[4]), chromosomes = chromosomes, cex = 1)
rainmut3si = plot_rainfall(vcfs[[6]], title = names(vcfs[6]), chromosomes = chromosomes, cex = 1)
rainwt1liver = plot_rainfall(vcfs[[7]], title = names(vcfs[7]), chromosomes = chromosomes, cex = 1)
rainwt2liver = plot_rainfall(vcfs[[8]], title = names(vcfs[8]), chromosomes = chromosomes, cex = 1)
rainwt3liver = plot_rainfall(vcfs[[10]], title = names(vcfs[10]), chromosomes = chromosomes, cex = 1)
rainwt2si = plot_rainfall(vcfs[[9]], title = names(vcfs[9]), chromosomes = chromosomes, cex = 1)
rainwt3si = plot_rainfall(vcfs[[11]], title = names(vcfs[11]), chromosomes = chromosomes, cex = 1)
blank = grid.rect(gp=gpar(col="white"))
rainfallplot = grid.arrange(rainwt1liver, rainmut1liver,blank,rainmut1si,
                            rainwt2liver,rainmut2liver, rainwt2si,rainmut2si, 
                            rainwt3liver,rainmut3liver,rainwt3si, rainmut3si, nrow=3, ncol =4)
# Save plot
#ggsave(paste(outdir,"extra/rainfallplot.pdf",sep=""), plot = rainfallplot, width = 20, height = 10)




# ---- TRANSCRIPTIONAL STRAND BIAS ----

#1 Get knowngenes table from UCSC for genome
genes_mm10 = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

#2 Make mutation count matrix with transcriptional information
mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_mm10)

#3 Strand bias
#3A Calculate
strand_counts = strand_occurrences(mut_mat_s, by=mousetype)
#3B Plot
strand_plot = plot_strand(strand_counts, mode = "relative")

#4 Log2 (transcribed/untranscribed)
#4A Calculate
strand_bias = strand_bias_test(strand_counts)
#4B Plot
strand_bias_plot = plot_strand_bias(strand_bias)

#5 Save plot
#ggsave(paste(outdir,"tcstrandbias.pdf",sep=""), plot = grid.arrange(strand_plot,strand_bias_plot, nrow=2, ncol =1), width = 10, height = 6)

#6 Statistical testing
#6A WT liver vs D- Liver (C>A (Transcribed(T)/total C>A;Untranscribed(U)/total C>A), C>G (T;U), C>T (T;U), T>A (T;U), T>C (T;U), T>G (T;U))
poisson.test(x = c(strand_counts[1,]$no_mutations,strand_counts[13,]$no_mutations), T = c(sum(strand_counts[1:2,]$no_mutations),sum(strand_counts[13:14,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[2,]$no_mutations,strand_counts[14,]$no_mutations), T = c(sum(strand_counts[1:2,]$no_mutations),sum(strand_counts[13:14,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[3,]$no_mutations,strand_counts[15,]$no_mutations), T = c(sum(strand_counts[3:4,]$no_mutations),sum(strand_counts[15:16,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[4,]$no_mutations,strand_counts[16,]$no_mutations), T = c(sum(strand_counts[3:4,]$no_mutations),sum(strand_counts[15:16,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[5,]$no_mutations,strand_counts[17,]$no_mutations), T = c(sum(strand_counts[5:6,]$no_mutations),sum(strand_counts[17:18,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[6,]$no_mutations,strand_counts[18,]$no_mutations), T = c(sum(strand_counts[5:6,]$no_mutations),sum(strand_counts[17:18,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[7,]$no_mutations,strand_counts[19,]$no_mutations), T = c(sum(strand_counts[7:8,]$no_mutations),sum(strand_counts[19:20,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[8,]$no_mutations,strand_counts[20,]$no_mutations), T = c(sum(strand_counts[7:8,]$no_mutations),sum(strand_counts[19:20,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[9,]$no_mutations,strand_counts[21,]$no_mutations), T = c(sum(strand_counts[9:10,]$no_mutations),sum(strand_counts[21:22,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[10,]$no_mutations,strand_counts[22,]$no_mutations), T = c(sum(strand_counts[9:10,]$no_mutations),sum(strand_counts[21:22,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[11,]$no_mutations,strand_counts[23,]$no_mutations), T = c(sum(strand_counts[11:12,]$no_mutations),sum(strand_counts[23:24,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[12,]$no_mutations,strand_counts[24,]$no_mutations), T = c(sum(strand_counts[11:12,]$no_mutations),sum(strand_counts[23:24,]$no_mutations))) #NS
#6B WT SI vs D- SI (C>A (Transcribed(T)/total C>A;Untranscribed(U)/total C>A), C>G (T;U), C>T (T;U), T>A (T;U), T>C (T;U), T>G (T;U))
poisson.test(x = c(strand_counts[25,]$no_mutations,strand_counts[37,]$no_mutations), T = c(sum(strand_counts[25:26,]$no_mutations),sum(strand_counts[37:38,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[26,]$no_mutations,strand_counts[38,]$no_mutations), T = c(sum(strand_counts[25:26,]$no_mutations),sum(strand_counts[37:38,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[27,]$no_mutations,strand_counts[39,]$no_mutations), T = c(sum(strand_counts[27:28,]$no_mutations),sum(strand_counts[39:40,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[28,]$no_mutations,strand_counts[40,]$no_mutations), T = c(sum(strand_counts[27:28,]$no_mutations),sum(strand_counts[39:40,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[29,]$no_mutations,strand_counts[41,]$no_mutations), T = c(sum(strand_counts[29:30,]$no_mutations),sum(strand_counts[41:42,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[30,]$no_mutations,strand_counts[42,]$no_mutations), T = c(sum(strand_counts[29:30,]$no_mutations),sum(strand_counts[41:42,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[31,]$no_mutations,strand_counts[43,]$no_mutations), T = c(sum(strand_counts[31:32,]$no_mutations),sum(strand_counts[43:44,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[32,]$no_mutations,strand_counts[44,]$no_mutations), T = c(sum(strand_counts[31:32,]$no_mutations),sum(strand_counts[43:44,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[33,]$no_mutations,strand_counts[45,]$no_mutations), T = c(sum(strand_counts[33:34,]$no_mutations),sum(strand_counts[45:46,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[34,]$no_mutations,strand_counts[46,]$no_mutations), T = c(sum(strand_counts[33:34,]$no_mutations),sum(strand_counts[45:46,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[35,]$no_mutations,strand_counts[47,]$no_mutations), T = c(sum(strand_counts[35:36,]$no_mutations),sum(strand_counts[47:48,]$no_mutations))) #NS
poisson.test(x = c(strand_counts[36,]$no_mutations,strand_counts[48,]$no_mutations), T = c(sum(strand_counts[35:36,]$no_mutations),sum(strand_counts[47:48,]$no_mutations))) #NS




# ---- GENOMIC DISTRIBUTION ----

#1 Choose mart
listEnsembl()
listMarts()
mart="ensembl"

#2 Choose dataset
# List datasets available from ensembl (for hg19 = GrCh37)
listDatasets(useEnsembl(biomart="regulation", GRCh = NULL))
# Multicell regulatory features
regulation_regulatory = useEnsembl(biomart="regulation", dataset="mmusculus_regulatory_feature", GRCh = NULL)
# list all possible filters
listFilters(regulation_regulatory)
# list all posible output attributes
listAttributes(regulation_regulatory)
# list all filter options for a specific attribute
filterOptions("regulatory_feature_type_name", regulation_regulatory)

#3 Get regions
#3A Promoters
promoter = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                 filters = "regulatory_feature_type_name", 
                 values = "Promoter", 
                 mart = regulation_regulatory)
promoter_g = reduce(GRanges(promoter$chromosome_name, IRanges(promoter$chromosome_start, promoter$chromosome_end)))
#3B Promoter flanking regions
promoter_flanking = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                          filters = "regulatory_feature_type_name", 
                          values = "Promoter Flanking Region", 
                          mart = regulation_regulatory)
promoter_flanking_g = reduce(GRanges(promoter_flanking$chromosome_name, IRanges(promoter_flanking$chromosome_start, promoter_flanking$chromosome_end))) 
#3C Enhancer
enhancer = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                 filters = "regulatory_feature_type_name", 
                 values = "Enhancer", 
                 mart = regulation_regulatory)
enhancer_g = reduce(GRanges(enhancer$chromosome_name, IRanges(enhancer$chromosome_start, enhancer$chromosome_end)))
#3D Genes
regionsg = list(genes_mm10)
names(regionsg) = c("Genes")
reduced_regionsg <- GenomicRanges::reduce(regionsg$Genes, ignore.strand=TRUE)

#4 Combine regions
regions = list(reduced_regionsg,promoter_g, promoter_flanking_g, enhancer_g)
names(regions) = c("Genes","Promoter", "Promoter flanking", "Enhancer")
regions = lapply(regions, function(x) rename_chrom(x))

#5 Get bedfiles (callable regions per sample)
surveyed_file = list.files(paste(indir,"callable_genome/bed/",sep=""), full.names = T)
surveyed_sample_names = c("Ercc1(-/D)3 Liver", "Ercc1(-/D)3 Small intestine", 
                          "Ercc1(-/D)2 Liver", "Ercc1(-/D)2 Small intestine", 
                          "Ercc1(-/D)1 Liver", "Ercc1(-/D)1 Small intestine",  
                          "WT1 Liver", 
                          "WT2 Liver","WT2 Small intestine",
                          "WT3 Liver", "WT3 Small intestine")
surveyed_list = bed_to_granges(surveyed_file, surveyed_sample_names)

#6 Calculate Log2(observed/expected) for all genomic regions
distr = genomic_distribution(vcfs, surveyed_list, regions)
distr_test = enrichment_depletion_test(distr, by = mousetype)

#7 Plot genomic distribution
genomic_dist_plot <- plot_enrichment_depletion(distr_test)
#ggsave(paste(outdir,"Genomicdistribution.pdf",sep=""), plot = genomic_dist_plot, width = 10, height = 6)

#8 Statistical testing
#8A Genes
aa_g <- poisson.test(x = c(distr_test[2,]$observed,distr_test[1,]$observed), T = c(distr_test[2,]$n_muts,distr_test[1,]$n_muts))$p.value #NS
cc_g <- poisson.test(x = c(distr_test[4,]$observed,distr_test[3,]$observed), T = c(distr_test[4,]$n_muts,distr_test[3,]$n_muts))$p.value #NS
#8B Promotor
aa_p <- poisson.test(x = c(distr_test[6,]$observed,distr_test[5,]$observed), T = c(distr_test[6,]$n_muts,distr_test[5,]$n_muts))$p.value #NS
cc_p <- poisson.test(x = c(distr_test[8,]$observed,distr_test[7,]$observed), T = c(distr_test[8,]$n_muts,distr_test[7,]$n_muts))$p.value #NS
#8C Promotor-flanking
aa_pf <- poisson.test(x = c(distr_test[10,]$observed,distr_test[9,]$observed), T = c(distr_test[10,]$n_muts,distr_test[9,]$n_muts))$p.value #NS
cc_pf <- poisson.test(x = c(distr_test[12,]$observed,distr_test[11,]$observed), T = c(distr_test[12,]$n_muts,distr_test[11,]$n_muts))$p.value #NS
#8D Enhancer
aa_e <- poisson.test(x = c(distr_test[14,]$observed,distr_test[13,]$observed), T = c(distr_test[14,]$n_muts,distr_test[13,]$n_muts))$p.value #NS
cc_e <- poisson.test(x = c(distr_test[16,]$observed,distr_test[15,]$observed), T = c(distr_test[16,]$n_muts,distr_test[15,]$n_muts))$p.value #NS
#8E Adjust for multiple testing
p.adjust(c(aa_g,cc_g), method = "fdr") #1.0000000 0.4340255
p.adjust(c(aa_p,cc_p), method = "fdr") #0.4977354 0.4977354
p.adjust(c(aa_pf,cc_pf), method = "fdr") #0.8849553 0.8849553
p.adjust(c(aa_e,cc_e), method = "fdr") #0.7417745 0.7417745
remove(aa_g,cc_g,aa_p,cc_p,aa_pf,cc_pf,aa_e,cc_e)
