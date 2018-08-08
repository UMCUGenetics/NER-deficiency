# @Date: 9 Feb 2017
# @Author: Francis Blokzijl
# @Description: Differential expression analysis ERCC1 delta min mouse organoids
# @Type: MUT & WT
# @Tissue: liver & small intestine
# @Raw data location:
# /hpc/cog_bioinf/cuppen/processed_data/nextseq/170117_NS500413_0253_AHC5NTBGX2_MOUSE

#biocLite("DESeq")
library(DESeq)
#library(DESeq2)
library(ggplot2)
#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
#install.packages("reshape")
library(reshape)
#install.packages("gplots")
library(gplots)
library(MutationalPatterns)
library(BSgenome)
# available.genomes()
ref_genome = "BSgenome.Mmusculus.UCSC.mm10" 
# biocLite(ref_genome)
library(ref_genome, character.only = TRUE)

# ---- RNASeq data ----

outdir = "~/Nextcloud/ERCC1/Results/RNAseq/DE/"

raw_counts_file = "~/Nextcloud/ERCC1/Data/RNAseq/170117_NS500413_0253_AHC5NTBGX2_MOUSE_readCounts_raw.txt"
RPKM_file = "~/Nextcloud/ERCC1/Data/RNAseq/170117_NS500413_0253_AHC5NTBGX2_MOUSE_readCounts_RPKM.txt"

countsTable = read.delim(raw_counts_file, header = T)
RPKM = read.delim(RPKM_file, header = T)

# add genes as row.names
rownames(countsTable) = countsTable[,1]
countsTable = countsTable[,-1]
rownames(RPKM) = RPKM[,1]
RPKM = RPKM[,-1]

# meta data
samples = colnames(countsTable)
conds = factor(c(rep("MUT",6), rep("WT",7)), levels = c("WT", "MUT"))
tissue = factor(c("LIVER", "SI", "LIVER", "SI", "LIVER", "SI", "LIVER","LIVER", "SI", "SI", "LIVER", "SI", "SI"))
cat = paste(conds, tissue, sep = "_")



# -------- COUNT NORMALIZATION ------------------
# DESeq2 function names are diferent? -> check vignette

# counts normalization
cds = newCountDataSet(countsTable, conds)
cds = estimateSizeFactors(cds)
normalized_counts = counts(cds, normalized=TRUE)
# Estimate variance
cds = estimateDispersions(cds)

# plot count histogram
df = melt(log2(normalized_counts))
colnames(df) = c("gene", 'sample', "log2normcounts")
plot_counts = ggplot(df, aes(x = log2normcounts)) +
  geom_histogram(binwidth=0.2, colour="black", fill="lightblue")

pdf(paste(outdir, "hist_nonzero_norm_counts.pdf", sep=""), 9, 6)
plot_counts
dev.off()

# plot RPKM histogram
df = melt(log2(RPKM))
colnames(df) = c("sample","log2RPKM")
plot_RPKM = ggplot(df, aes(x = log2RPKM)) +
  geom_histogram(binwidth=0.2, colour="black", fill="lightblue")

pdf(paste(outdir, "hist_nonzero_RPKM.pdf", sep=""), 9, 6)
plot_RPKM
dev.off()


# ------- SAMPLE CORRELATION ---------
# discard genes with < 10 counts over all samples
expressed = which(rowSums(normalized_counts) > 10)

# spearman correlation between samples based on expressed genes
cor_samples = cor(normalized_counts[expressed,], method = 'spearman')
plot_cor = heatmap(cor_samples, scale = 'none', col = rev(heat.colors(256)))

pdf(paste(outdir, "correlation_samples.pdf", sep=""), 11, 11)
heatmap(cor_samples, scale = 'none', col = rev(heat.colors(256)))
dev.off()


# ------ DIFFERENTIAL EXPRESSION ------ 

# DE WT VS MUTANT OVERALL (both liver & small intestine)

# Test for differential expression
# H0 : the measurements come from the same distribution and the gene is being expressed 
# at the same level across conditions

# general function for differential expression calculation
# annotate with gene symbols & description
diff_expr = function(counts, conds, map)
{
  samples = colnames(counts)
  # counts normalization
  cds = newCountDataSet(counts, conds)
  cds = estimateSizeFactors(cds)
  normalized_counts = counts(cds, normalized=TRUE)
  # Estimate variance
  cds = estimateDispersions(cds)
  DE = nbinomTest(cds, levels(conds)[1], levels(conds)[2])
  # add normalized counts
  DE = merge(normalized_counts, DE, by.x= "row.names", by.y = "id")
  # recalculate log2foldchange with pseudocount to avoid Inf
  DE$log2FoldChange_pseudo= log2((DE$baseMeanB + 1) / (DE$baseMeanA+1))
  # find gene symbol and description
  cols = c("SYMBOL","GENENAME")
  DE_genes = select(map, keys=DE$Row.names, columns=cols, keytype="ENSEMBL")
  DE = merge(DE, DE_genes, by.x = "Row.names", by.y = "ENSEMBL")
  # remove duplicate genes (which arises during gene conversion)
  DE = DE[which(!duplicated(DE$Row.names)),]
  # sort result on adjusted pvalue
  DE = DE[order(DE$padj),]
  return(DE)
}

# DE MUT & WT
DE_overall = diff_expr(countsTable, conds, org.Mm.eg.db)
# DE MUT & WT liver
liver_samples = samples[tissue == "LIVER"]
liver_conds = conds[tissue == "LIVER"]
DE_liver = diff_expr(countsTable[,liver_samples], liver_conds, org.Mm.eg.db)
# DE MUT & WT SI
SI_samples = samples[tissue == "SI"]
SI_conds = conds[tissue == "SI"]
DE_SI = diff_expr(countsTable[,SI_samples], SI_conds, org.Mm.eg.db)
# DE liver & SI
DE_tissues = diff_expr(countsTable, tissue, org.Mm.eg.db)

#DE WT
wt_samples = samples[conds == "WT"]
wt_conds = tissue[conds == "WT"]
DE_WT = diff_expr(countsTable[,wt_samples],wt_conds,org.Mm.eg.db)


# Find interesting DE genes
threshold = 2
DE_overall_int = subset(DE_overall, DE_overall$padj < 0.05 & abs(DE_overall$log2FoldChange_pseudo) > threshold)
DE_tissues_int = subset(DE_tissues, DE_tissues$padj < 0.05 & abs(DE_tissues$log2FoldChange_pseudo) > threshold)
DE_liver_int = subset(DE_liver, DE_liver$padj < 0.05 & abs(DE_liver$log2FoldChange_pseudo) > threshold)
DE_SI_int = subset(DE_SI, DE_SI$padj < 0.05 & abs(DE_SI$log2FoldChange_pseudo) > threshold)
DE_WT_int = subset(DE_WT, DE_WT$padj < 0.05 & abs(DE_WT$log2FoldChange_pseudo) > threshold)

# write output
write.table(DE_liver_int, file = paste(outdir, "DE_liver_interesting.txt", sep=""), sep = "\t", quote = F, row.names = F)
write.table(DE_liver, file = paste(outdir, "DE_liver_all.txt", sep =""), sep = "\t", quote = F, row.names = F)
write.table(DE_SI_int, file = paste(outdir, "DE_SI_interesting.txt", sep=""), sep = "\t", quote = F, row.names = F)
write.table(DE_SI, file = paste(outdir, "DE_SI_all.txt", sep =""), sep = "\t", quote = F, row.names = F)


# ------- PLOT EXPRESSION -------

# Function to plot expression for a specific gene
plot_gene_expression = function(gene_symbol)
{
  df = subset(DE_overall, DE_overall$SYMBOL == gene_symbol)
  df = melt(df[,2:14])
  df$Tissue = tissue
  df$Type = conds
  colnames(df)[1] = "sample"
  colnames(df)[2] = "Normalised_counts"
  
  plot = ggplot(df, aes(x= Type, y= Normalised_counts, fill=Tissue)) +
    geom_boxplot() +
    geom_point() +
    facet_grid(. ~ Tissue) +
    theme(legend.position="none") +
    ggtitle(gene_symbol)
  return(plot)
}

plot_gene_expression_WT = function(gene_symbol)
{
  df = subset(DE_WT, DE_WT$SYMBOL == gene_symbol)
  df = melt(df[,2:8])
  df$Tissue = wt_conds
#  df$Type = conds
  colnames(df)[1] = "sample"
  colnames(df)[2] = "Normalised_counts"
  
  plot = ggplot(df, aes(x= Tissue, y= Normalised_counts, fill=Tissue)) +
    geom_boxplot() +
    geom_point() +
 #   facet_grid(. ~ Tissue) +
    theme(legend.position="none") +
    expand_limits(y=0) +
    ggtitle(gene_symbol)
  return(plot)
}




# Plot ERCC1, XPF & Lgr5 expression
pdf(paste(outdir, "ERCC1_expression.pdf", sep=""), 6, 4)
plot_gene_expression("Ercc1")
dev.off()

pdf(paste(outdir, "XPF_expression.pdf", sep=""), 6, 4)
plot_gene_expression("Ercc4")
dev.off()

pdf(paste(outdir, "lgr5_expression.pdf", sep=""), 6, 4)
plot_gene_expression("Lgr5")
dev.off()

# interesting repair genes gene list
gene_list = c("Xpc", "Rad23b", "Ercc8", "Ercc6", "Xpa", "Ercc3", "Ercc2", 
              "Ercc1", "Ercc4", "Ercc5", "Xrcc1", "Xrcc2", "Xrcc3", "Xrcc4",
              "Xrcc5", "Xrcc6", "Brca1", "Brca2", "Trp53")
# oxidative stress gene list
gene_list_oxstress = c("Ogg1", "Mutyh", "Nudt1")
# apoptosis gene list
gene_list_apoptosis = c("Bik", "Bax", "Bad", "Fas", "Fis1", "Fxn", "Gpx1", "Igf1", "Pmaip1", "Tnfsf10", "Casp8")

                       
# plot expression of gene lists
pdf(paste(outdir, "repair_gene_list_expression.pdf", sep=""), 6, 4)
for(i in 1:length(gene_list)){
  print(plot_gene_expression(gene_list[i]))
}
dev.off()

pdf(paste(outdir, "apoptosis_gene_list_expression.pdf", sep=""), 6, 4)
for(i in 1:length(gene_list_apoptosis)){
  print(plot_gene_expression(gene_list_apoptosis[i]))
}
dev.off()

pdf(paste(outdir, "oxstress_gene_list_expression.pdf", sep=""), 6, 4)
for(i in 1:length(gene_list_oxstress)){
  print(plot_gene_expression(gene_list_oxstress[i]))
}
dev.off()

#NU MEER DNA REPAIR GERICHT
dnarepair <- read.table("~/Nextcloud/ERCC1/Data/RNAseq/ercc1_dnarepair.txt", header = F)
dnarepair <- as.character(dnarepair$V1)
dnarepair <- toupper(dnarepair)
dnarepair <- dnarepair[-which(dnarepair == "ERCC1")]
ner_repair_lower <- c("Xpc","Ercc8", "Ercc6", "Xpa", "Ercc3", "Ercc2", 
                      "Ercc1", "Ercc4", "Ercc5")

outdir <- "~/Nextcloud/ERCC1/Results/RNAseq/"
pdf(paste(outdir, "ner_expression.pdf", sep=""), 6, 4)
for(i in 1:length(ner_repair_lower)){
  print(plot_gene_expression_WT(ner_repair_lower[i]))
}
dev.off()


#DE DNA repair

DE_WT_NER <- DE_WT[1,]
DE_WT_NER <- DE_WT_NER[-1,]
DE_WT_REP <- DE_WT_NER

DE_SI_REP_ALL<- DE_WT_NER
DE_LIVER_REP_ALL<- DE_WT_NER
DE_REP_ALL <- DE_WT_NER

p.adjust(c(DE_SI[which(DE_SI$SYMBOL == "Ercc1"),]$pval,DE_liver[which(DE_liver$SYMBOL == "Ercc1"),]$pval), method = "fdr")

for(i in ner_repair_lower) {
  temp <- subset(DE_WT,DE_WT$SYMBOL == i)
  DE_WT_NER <- rbind(DE_WT_NER,temp)
}
DE_WT_NER$padj <- p.adjust(c(DE_WT_NER$pval), method = "fdr")

for(i in dnarepair) {
  temp <-subset(DE_SI,toupper(DE_SI$SYMBOL) == i)
  DE_SI_REP_ALL <- rbind(DE_SI_REP_ALL,temp)
}
DE_SI_REP_ALL$padj <- p.adjust(c(DE_SI_REP_ALL$pval), method = "fdr")

for(i in dnarepair) {
  temp <-subset(DE_liver,toupper(DE_liver$SYMBOL) == i)
  DE_LIVER_REP_ALL <- rbind(DE_LIVER_REP_ALL,temp)
}
DE_LIVER_REP_ALL$padj <- p.adjust(c(DE_LIVER_REP_ALL$pval), method = "fdr")

for(i in dnarepair) {
  temp <-subset(DE_overall,toupper(DE_overall$SYMBOL) == i)
  DE_REP_ALL <- rbind(DE_REP_ALL,temp)
}
DE_REP_ALL$padj <- p.adjust(c(DE_REP_ALL$pval), method = "fdr")

output_NER <- DE_WT_NER[,c(1,17,10,11,13,16,14,15)]
colnames(output_NER) <- c("Ensembl gene ID","Gene symbol","Mean expression (Liver)","Mean expression (SI)","log2FoldChange","log2FoldChange_pseudo","p-value","Adjusted p-value")
rownames(output_NER) <- c()
write.table(file = "~/Nextcloud/ERCC1/Results/RNAseq/ner.txt", x = output_NER, sep = "\t",row.names = FALSE)

output_REP_SI <-DE_SI_REP_ALL[,c(1,17,10,11,13,16,14,15)]
colnames(output_REP_SI) <- c("Ensembl gene ID","Gene symbol","Mean expression (WT SI)","Mean expression (Ercc1D/- SI)","log2FoldChange (SI)","log2FoldChange_pseudo (SI)","p-value (SI)","Adjusted p-value (SI)")
rownames(output_REP_SI) <- c()
write.table(file = "~/Nextcloud/ERCC1/Results/RNAseq/repair_si.txt", x=output_REP_SI, sep = "\t",row.names = FALSE)

output_REP_Liver <- DE_LIVER_REP_ALL[,c(1,16,9,10,12,15,13,14)]
colnames(output_REP_Liver) <- c("Ensembl gene ID","Gene symbol","Mean expression (WT Liver)","Mean expression (Ercc1D/- Liver)","log2FoldChange (Liver)","log2FoldChange_pseudo (Liver)","p-value (Liver)","Adjusted p-value (Liver)")
rownames(output_REP_Liver) <- c()
write.table(file = "~/Nextcloud/ERCC1/Results/RNAseq/repair_liver.txt", x=output_REP_Liver, sep = "\t",row.names = FALSE)

fraction_NER <- nrow(DE_WT_NER[which(DE_WT_NER$padj < 0.05 & DE_WT_NER$log2FoldChange_pseudo < 0),])/nrow(DE_WT_NER)
fraction_NER_diff <- nrow(DE_WT_NER[which(DE_WT_NER$padj < 0.05),])/nrow(DE_WT_NER)
fraction_REP <- nrow(DE_WT_REP[which(DE_WT_REP$padj < 0.05 & DE_WT_REP$log2FoldChange_pseudo < 0),])/nrow(DE_WT_REP)
fraction_REP_diff <- nrow(DE_WT_REP[which(DE_WT_REP$padj < 0.05),])/nrow(DE_WT_REP)
fraction_total <- nrow(DE_WT[which(DE_WT$padj < 0.05 & DE_WT$log2FoldChange_pseudo < 0),])/nrow(DE_WT)
fraction_total_diff <- nrow(DE_WT[which(DE_WT$padj < 0.05),])/nrow(DE_WT)

pval_NER <- poisson.test(x=c(nrow(DE_WT_NER[which(DE_WT_NER$padj < 0.05 & DE_WT_NER$log2FoldChange_pseudo < 0),]),nrow(DE_WT_NER)),
             T=c(nrow(DE_WT[which(DE_WT$padj < 0.05 & DE_WT$log2FoldChange_pseudo < 0),]),nrow(DE_WT)))$p.value
pval_REP <- poisson.test(x=c(nrow(DE_WT_REP[which(DE_WT_REP$padj < 0.05 & DE_WT_REP$log2FoldChange_pseudo < 0),]),nrow(DE_WT_REP)),
                         T=c(nrow(DE_WT[which(DE_WT$padj < 0.05 & DE_WT$log2FoldChange_pseudo < 0),]),nrow(DE_WT)))$p.value
p.adjust(c(pval_NER,pval_REP),method = "fdr") #0.02729068 0.02729068



#get epxressed genes

nrow(countsTable[which(countsTable[,14] >1),])

expressed_genes <- c(14745,16353,16447,16190,16027,16866,15977,15713,15862,16096,16028,16383,16277)
mean(expressed_genes) #16074
sd(expressed_genes) #496
