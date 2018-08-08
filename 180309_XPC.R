# @Date: 30 October 2017
# @Modified: 08 March 2018
# @Author: Myrthe Jager
# @Description: Mutational patterns in XPC WT and XPC KO organoids
# Abbreviations: N = normal




# ---- GET STARTED ----

#1 Install & load required packages
library(MutationalPatterns)
library(BSgenome)
library(reshape2)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("gridExtra")
library(grid)
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
library(gridExtra)

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/XPC/"

#3 Functions
## Function 1: Calculate intermutation distance ##
mut_dist = function(vcf) {
  # mutation characteristics
  type = loc = dist = chrom = previous = c()
  
  # for each chromosome
  for(i in 1:length(chromosomes))
  {
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    n = length(chr_subset)
    if(n<=1){next}
    type = c(type, mut_type(chr_subset)[-1])
    loc = c(loc, (start(chr_subset))[-1])# + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    chrom = c(chrom, rep(chromosomes[i],n-1))
  }
  
  data = data.frame(type = type,
                    location = loc,
                    distance = dist,
                    chromosome = chrom)
  return(data)
}

## Function 2: Rename CHROM to UCSC style ##
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

## Function 3: Convert bed to Granges ##
bed_to_granges = function(bed_files, region_names) { 
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

#4 Install and load human reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 metadata
sample_name <- c("N1","N2","N3","XPC KO")
type_iv_xpc <- c("normal","normal","normal","XPC KO")




# ---- LOAD DATA: SNV VCF ----

#1 List vcf files: filtered SNVs
vcf_file = list.files("~/surfdrive/shared/ERCC1/Data/XPC/SNV/", full.names = TRUE)

#2 read vcfs as granges
vcf = read_vcfs_as_granges(vcf_file, sample_name, genome = ref_genome)

#3 Check if chromosome names are uniform
all(seqlevels(vcf[[1]]) %in% seqlevels(get(ref_genome)))

#4 only select autosomal chromosomes
auto = extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")
vcf = lapply(vcf, function(x) keepSeqlevels(x, auto))




# ---- FIND DOUBLE BASE SUBSTITUTIONS ----

#1 get cumulative sum of chromosome lengths
chr_length = seqlengths(vcf[[1]])
chromosomes = names(chr_length)
chr_cum = c(0, cumsum(as.numeric(chr_length)))

#2 Calculate intermutation distance
mut_dist_vcfs = lapply(vcf, function(x) mut_dist(x))

#3 Generate table to check double nucleotides manually
#3A Extract double nucleotides
confirm <- mut_dist_vcfs$N1[which(mut_dist_vcfs$N1[,3] == 1),]
confirm <- rbind(confirm,mut_dist_vcfs$N2[which(mut_dist_vcfs$N2[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$N3[which(mut_dist_vcfs$N3[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`XPC KO`[which(mut_dist_vcfs$`XPC KO`[,3] == 1),])
#3B Put in a table and add sample type
dinucleotide <- data.frame(
  chr = confirm$chromosome,
  pos = confirm$location-1,
  sample <- c(rep("N1",nrow(mut_dist_vcfs$N1[which(mut_dist_vcfs$N1[,3] == 1),])),
              rep("N2",nrow(mut_dist_vcfs$N2[which(mut_dist_vcfs$N2[,3] == 1),])),
              rep("N3",nrow(mut_dist_vcfs$N3[which(mut_dist_vcfs$N3[,3] == 1),])),
              rep("XPC KO",nrow(mut_dist_vcfs$`XPC KO`[which(mut_dist_vcfs$`XPC KO`[,3] == 1),]))
  )
)
colnames(dinucleotide) <- c("chr","pos","sample")
remove(confirm)
#3C Write table
#write.table(dinucleotide,paste(outdir,"Dinucleotides_XPC_FPandTP.txt",sep = ""),sep="\t", row.names = FALSE, col.names = TRUE)
#3D These mutations were checked manually in IGV




# ---- SNV LOAD: SINGLE AND DOUBLE ----

#1 Import table: Callable genome
xpccorrectedmutnumber <- read.delim2(paste(indir,"XPC/Callable_XPC.txt", sep = ""))
xpccorrectedmutnumber$type <- type_iv_xpc
xpccorrectedmutnumber$sample <- sample_name
xpccorrectedmutnumber$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(xpccorrectedmutnumber$Surveyed.percentage)))

#2 Add nucleotide substitutions counts (single + double!)
xpccorrectedmutnumber$all.uncorr.SNV <- NA
for (i in vcf_file) {
  df.temp <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
  snv_number <- nrow(df.temp)
  
  # Add to callable genome table
  xpccorrectedmutnumber[which(toupper(xpccorrectedmutnumber$Sample) == toupper(sample)),]$all.uncorr.SNV <- snv_number
}
remove(df.temp,i,sample,snv_number)

#3 Double nucleotide substitution counts
#3A Import table: IGV-verfied double base substitutions
checkeddouble <- read.delim2(paste(indir,"XPC/xpc_dinucleotides_igv.txt",sep=""))
checkeddouble <- checkeddouble[which(checkeddouble$true == "TRUE"),]
#3B Counts double nucleotide substitutions and add counts to table
xpccorrectedmutnumber$uncorr.double <- c(as.numeric(length(which(checkeddouble$sample == "N1"))),
                                            as.numeric(length(which(checkeddouble$sample == "N2"))),
                                            as.numeric(length(which(checkeddouble$sample == "N3"))),
                                            as.numeric(length(which(checkeddouble$sample == "XPC KO"))))

#4 Correct SNVs for double nucleotide substitutions
xpccorrectedmutnumber$uncorr.snv <- xpccorrectedmutnumber$all.uncorr.SNV-(2*xpccorrectedmutnumber$uncorr.double)

#5 Calculate the corrected number of nucleotide substitutions (extrapolate to autosomal genome)
#5A Single nucleotide substitutions
xpccorrectedmutnumber$corr.snv <- xpccorrectedmutnumber$uncorr.snv/xpccorrectedmutnumber$Surveyed.percentage*100
#5B Double nucleotide substitutions
xpccorrectedmutnumber$corr.double <- xpccorrectedmutnumber$uncorr.double/xpccorrectedmutnumber$Surveyed.percentage*100

#6 Weeks in culture
xpccorrectedmutnumber$weeks <- c(rep(20.57143,3),10.28571)

#7 Calculate the number of double nucleotide substitutions per week
#7A Single nucleotide substitutions
xpccorrectedmutnumber$corr.snv.perweek <- xpccorrectedmutnumber$corr.snv/xpccorrectedmutnumber$weeks
#7B Double nucleotide substitutions
xpccorrectedmutnumber$corr.double.perweek <- xpccorrectedmutnumber$corr.double/xpccorrectedmutnumber$weeks
#7C Write table with mutation numbers
#write.table(xpccorrectedmutnumber,file = paste(outdir,"xpc_mutationnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)

#8 Calculate increase per week
#8A SNVs
xpccorrectedmutnumber[which(xpccorrectedmutnumber$Sample == "STE0076XPCSC44"),]$corr.snv.perweek/mean(xpccorrectedmutnumber[which(xpccorrectedmutnumber$Sample != "STE0076XPCSC44"),]$corr.snv.perweek)
#8A Double nucleotide substitutions
xpccorrectedmutnumber[which(xpccorrectedmutnumber$Sample == "STE0076XPCSC44"),]$corr.double.perweek/mean(xpccorrectedmutnumber[which(xpccorrectedmutnumber$Sample != "STE0076XPCSC44"),]$corr.double.perweek)

#9 Plot
#9A Single nucleotide substitutions
xpc_snv_numberplot <- ggplot(data=xpccorrectedmutnumber, aes(x=sample, y=corr.snv.perweek, fill= type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab("Xpc") + 
  ylab("Base substitutions\nper autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=100) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())
#ggsave(paste(outdir,"xpc_SNVnumber.pdf",sep=""), xpc_snv_numberplot, width = 8, height = 8)
#9B Double nucleotide substitutions
xpc_double_numberplot <- ggplot(data=xpccorrectedmutnumber, aes(x=sample, y=corr.double.perweek, fill= type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab("Xpc") + 
  ylab("Double base substitutions\nper autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=3) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())
#ggsave(paste(outdir,"xpc_DOUBLEnumber.pdf",sep=""), xpc_double_numberplot, width = 8, height = 8)




# ---- MUTATION SPECTRA ----

#1 Mutation spectra (6 types)
type_occurrence <- mut_type_occurrences(vcf, ref_genome)

#2 Plot mutation spectra
plot_spectrum(type_occurrence, by = type_iv_xpc,legend = T)

#3 Mutational profiles (96 types)
mut_matrix_xpc = mut_matrix(vcf_list = vcf, ref_genome = ref_genome)




# ---- REFIT SAMPLES & SIGNATURES ----

#1 Signatures for refitting
#1A Download signatures from pan-cancer study Alexandrov et al.
cancer_signatures = read.table(paste(indir,"SNV/cancersignatures.txt",sep=""), sep="\t")
#1B Reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
#1C Only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
#1D ERCC1 sigs
ercc1_signatures = read.table(file = paste(indir,"SNV/NMF_2signatures.txt",sep=""), sep = "\t")
ercc1_signatures <- as.matrix(ercc1_signatures)

#2 Fit to signatures
#2A Cosmic cancer signatures
fit_res_xpc_cancersigs = fit_to_signatures(mut_matrix_xpc, cancer_signatures)
#2B Ercc1 signatures
fit_res_xpc_ercc1sigs <- fit_to_signatures(mut_matrix_xpc, ercc1_signatures)

#3 Contribution per week
#3A Get absolute contribution
perweek_cancersigs <- fit_res_xpc_cancersigs$contribution
perweek_ercc1 <- fit_res_xpc_ercc1sigs$contribution
#3B Add number of weeks
perweek_cancersigs <- rbind(perweek_cancersigs,xpccorrectedmutnumber$weeks)
perweek_ercc1 <- rbind(perweek_ercc1,xpccorrectedmutnumber$weeks)
#3C Calculate number of mutations per signature per week
for (i in 1:ncol(perweek_cancersigs)){
  for(j in 1:nrow(perweek_cancersigs)) {
    perweek_cancersigs[j,i] <- perweek_cancersigs[j,i]/perweek_cancersigs[31,i]
  }
}
perweek_cancersigs<-perweek_cancersigs[1:30,]
for (i in 1:ncol(perweek_ercc1)){
  for(j in 1:nrow(perweek_ercc1)) {
    perweek_ercc1[j,i] <- perweek_ercc1[j,i]/perweek_ercc1[3,i]
  }
}
perweek_ercc1 <- perweek_ercc1[1:2,]

#4 Plot contribution of all signatures per week
#4A Cosmic signatures
contribution.plot.cancersigs <- plot_contribution(perweek_cancersigs, coord_flip = F, signatures = cancer_signatures, mode = "absolute")
#ggsave(paste(outdir,"xpc_cancersigs_perweek.pdf",sep=""), plot = contribution.plot.cancersigs,width = 8, height = 8)
#4B Ercc1 signatures
contribution.plot.ercc1sigs <- plot_contribution(perweek_ercc1, coord_flip = F, signatures = ercc1_signatures, mode = "absolute")

#5 Contribution of signature 8 per week
#5A Generate data frame
sig8<- data.frame(sig8 = perweek_cancersigs[8,],
  name = sample_name,
  type = type_iv_xpc,
  surveyed = xpccorrectedmutnumber$Surveyed.percentage)
#5B Extrapolate to autosomal genome
sig8$corrected <- sig8$sig8/sig8$surveyed*100
#5C Plot
plot_sig8<- ggplot(data=sig8, aes(x=name, y=corrected, fill= type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab("Xpc") + 
  ylab("Signature 8 base substitutions\nper autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=40) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())
#ggsave(paste(outdir,"xpc_sig8_perweek.pdf",sep=""), plot = plot_sig8,width = 8, height = 8)

#6 Signature contribution per week: plot all cosmic signatures separately
#6A Extrapolate to autosomal genome
corrected_perweek_cancersigs <- data.frame()
for(i in 1:nrow(perweek_cancersigs)) {
  corrected <- perweek_cancersigs[i,]/xpccorrectedmutnumber$Surveyed.percentage*100
  corrected_perweek_cancersigs <- rbind(corrected_perweek_cancersigs,corrected)
  colnames(corrected_perweek_cancersigs) <- colnames(perweek_cancersigs)
  remove(corrected)
}
remove(i)
row.names(corrected_perweek_cancersigs) <- row.names(perweek_cancersigs)
#6B generate table with mean and STdev for the WT
sig.contribution <- data.frame()
for (i in 1:nrow(corrected_perweek_cancersigs)) {
  wt.mean <- sum(corrected_perweek_cancersigs[i,1:3])/3
  wt.stdev <- sd(corrected_perweek_cancersigs[i,1:3])
  df <- data.frame(xpc.wt= wt.mean,
                   wt.stdev = wt.stdev,
                   xpc.ko = corrected_perweek_cancersigs[i,4],
                   mut.stdev = NA)
  row.names(df) <- row.names(corrected_perweek_cancersigs)[i]
  sig.contribution <- rbind(sig.contribution,df)
  remove(wt.mean, wt.stdev,df)
}
remove(i)
#6C Plot 
sig.contribution$Signature <- factor(row.names(sig.contribution),levels = row.names(sig.contribution))
m.sig.contribution <- melt(sig.contribution[,c(1,3,5)],id.vars = 'Signature')
m.sig.contribution$stdev <- melt(sig.contribution[,c(2,4,5)],id.vars = 'Signature')$value
m.sig.contribution$max <- m.sig.contribution$value+m.sig.contribution$stdev 
m.sig.contribution$min <- m.sig.contribution$value-m.sig.contribution$stdev 
m.sig.contribution[which(m.sig.contribution$min < 0),]$min <- 0
#type.forplot = factor(x=c("XPC WT","XPC KO"),levels = c("XPC WT","XPC KO"))
plot.sig.all <-  ggplot(data = m.sig.contribution, aes(x = Signature, y = value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30")) +
  geom_errorbar(aes(ymin=min, ymax=max),
              width=.03,                    # Width of the error bars
              position=position_dodge(.7)) +
  ylim(0,40) +
  ylab("Base substitutions\nper autosomal genome per week\n") +
  theme(axis.text=element_text(size = rel(0.6)), 
        axis.title=element_text(size = rel(0.8)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        axis.line.x = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"))
#ggsave(paste(outdir,"all_sigs_separate.pdf",sep=""),plot = plot.sig.all,width = 6.5, height = 2) 

#7 Percentage signature 8
abs_si8 <- sig8[4,]$corrected-mean(sig8[1:3,]$corrected)
abs_total <- (xpccorrectedmutnumber$all.uncorr.SNV/xpccorrectedmutnumber$Surveyed.percentage*100)/xpccorrectedmutnumber$weeks
abs_total <- abs_total[4]-mean(abs_total[1:3])
abs_si8/abs_total*100 #39.85717%

#8 Contribution of ercc1 signatures per week
#8A Generate data frame
ercc1_pw<- data.frame(siga = perweek_ercc1[1,],
                      sigb = perweek_ercc1[2,],
                  name = sample_name,
                  type = type_iv_xpc,
                  surveyed = xpccorrectedmutnumber$Surveyed.percentage)
#8B Extrapolate to autosomal genome
ercc1_pw$corrected_siga <- ercc1_pw$siga/ercc1_pw$surveyed*100
ercc1_pw$corrected_sigb <- ercc1_pw$sigb/ercc1_pw$surveyed*100
#8C Generate data frame for plotting
ercc1_sig_mutnumber <- data.frame(
  name = c(as.character(ercc1_pw$name),as.character(ercc1_pw$name)),
  type =  c(as.character(ercc1_pw$type),as.character(ercc1_pw$type)),
  mutationtype = c(rep("Signature A",4),rep("Signature B",4)),
  mutations.perweek = c(ercc1_pw$corrected_siga,ercc1_pw$corrected_sigb)
)
#8D Plot
plot_ercc1_sigs<- ggplot(data=ercc1_sig_mutnumber, aes(x=name, y=mutations.perweek, fill= mutationtype)) +
  geom_bar(stat="identity", width = 0.7) +
  xlab("Xpc") + 
  ylab("Absolute contribution\n(base substitutions\nper autosomal genome\nper week)\n") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=100) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())
#ggsave(paste(outdir,"xpc_ercc1sigs_perweek.pdf",sep=""), plot = plot_ercc1_sigs,width = 8, height = 8)




# ---- STRAND BIAS ----

#1 Get knowngenes table from UCSC for genome
genes_hg19 = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

#2 Make mutation count matrix with transcriptional information
mut_mat_s_xpc = mut_matrix_stranded(vcf, ref_genome, genes_hg19)

#3 Strand bias
#3A Calculate
strand_counts_xpc = strand_occurrences(mut_mat_s_xpc, by = type_iv_xpc)
#3B Plot
strand_plot_xpc = plot_strand(strand_counts_xpc, mode = "relative")

#4 Log2 (transcribed/untranscribed)
#4A Calculate
strand_bias_xpc = strand_bias_test(strand_counts_xpc)
#4B Plot
strand_bias_plot_xpc = plot_strand_bias(strand_bias_xpc)

#5 Save plot
#ggsave(paste(outdir,"xpc_strandbias.pdf",sep=""), plot=grid.arrange(strand_plot_xpc,strand_bias_plot_xpc, nrow=2, ncol =1))




# ---- GENOMIC DISTRIBUTION ----

#1 Choose mart
listEnsembl()
listMarts()
mart="ensembl"

#2 Choose dataset
listDatasets(useEnsembl(biomart="regulation", GRCh = 37))
# Multicell regulatory features for hg19
regulation_regulatory = useEnsembl(biomart="regulation", dataset="hsapiens_regulatory_feature", GRCh = 37)
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
regionsg = list(genes_hg19)
names(regionsg) = c("Genes")
reduced_regionsg <- GenomicRanges::reduce(regionsg$Genes, ignore.strand=TRUE)

#4 Combine regions
regions_final_xpc = list(reduced_regionsg ,promoter_g, promoter_flanking_g, enhancer_g)
names(regions_final_xpc)= c("Genes","Promoter", "Promoter flanking", "Enhancer")
regions_final_xpc = lapply(regions_final_xpc,function(x) rename_chrom(x))

#5 Get bedfiles (callable regions per sample)
surveyed_file_xpc = list.files(paste(indir,"XPC/final_callableLoci/",sep=""), full.names = T)
surveyed_list_xpc = bed_to_granges(surveyed_file_xpc, sample_name)

#6 Calculate Log2(observed/expected) for all genomic regions
distr_final_xpc = genomic_distribution(vcf, surveyed_list_xpc, regions_final_xpc)
distr_test_xpc = enrichment_depletion_test(distr_final_xpc, by = type_iv_xpc)
distr_test_xpc_persample = enrichment_depletion_test(distr_final_xpc)

#7 Plot genomic distribution
genomic_dist_plot <- plot_enrichment_depletion(distr_test_xpc)
#ggsave(paste(outdir,"xpc_Genomicdistribution.pdf",sep=""), plot = genomic_dist_plot, width = 10, height = 6)