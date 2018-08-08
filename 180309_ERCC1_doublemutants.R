# @Date: 30 October 2017
# @Modified: 08 March 2018 and 09 March 2018
# @Author: Myrthe Jager & Francis Blokzijl
# @Description: double base substitution analysis ERCC1 delta min mouse organoids
# Abbreviations: MUT = Ercc1 delta min; WT = Wild-type, SI = Small intestine; SNV = single nucleotide variant (point mutation)




# ---- GET STARTED ----

#1 Install & load required packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome")
library(BSgenome)
#biocLite("MutationalPatterns")
library(MutationalPatterns)
#install.packages("ggplot2")
library(ggplot2)

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/SNV/Dinucleotide/"

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

## Function 2: Calculate intermutation distance ##
mut_dist = function(vcf){
  # mutation characteristics
  type = loc = dist = chrom = previous = c()
  
  # for each chromosome
  for(i in 1:length(chromosomes)){
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

## Function 3: Calculate whiskers based on sd ##
remove.negative.sd <- function(df) {
  df$pos <- df[,3] + df$sd
  df$neg <- df[,3] - df$sd
  if (length(df[which(df$neg < 0),]$neg) >0) {
    df[which(df$neg < 0),]$neg <- 0
  }
  return(df)
}

#4 Install and load mouse reference genome
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
#biocLite(ref_genome)
library(ref_genome, character.only = T)

#5 metadata
vcf_type = c(rep("MUT", 6), rep("WT", 5))
vcf_type = factor(vcf_type, levels = c("WT", "MUT"))
vcf_tissue = c(rep(c("Liver","Small intestine"),3),"Liver",rep(c("Liver","Small intestine"),2))
vcf_cond = paste(vcf_type, vcf_tissue, sep = "_")
sample_names = c("Ercc1(-/D)1 Liver", "Ercc1(-/D)1 Small intestine", 
                 "Ercc1(-/D)2 Liver", "Ercc1(-/D)2 Small intestine", 
                 "Ercc1(-/D)3 Liver", "Ercc1(-/D)3 Small intestine",  
                 "WT1 Liver", 
                 "WT2 Liver","WT2 Small intestine",
                 "WT3 Liver", "WT3 Small intestine")




# ---- LOAD DATA: SNV VCF ----

#1 List vcf files: filtered SNVs
vcf_files = list.files(paste(indir,"SNV/VCF/",sep=""), full.names = T)

#2 read vcfs as granges
vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)

#3 Check if chromosome names are uniform
all(seqlevels(vcfs[[1]]) %in% seqlevels(get(ref_genome)))

#4 only select autosomal chromosomes
auto = extractSeqlevelsByGroup(species="Mus_musculus", style="UCSC", group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))




# ---- FIND DOUBLE BASE SUBSTITUTIONS ----

#1 get cumulative sum of chromosome lengths
chr_length = seqlengths(vcfs[[1]])
chromosomes = names(chr_length)
chr_cum = c(0, cumsum(as.numeric(chr_length)))

#2 Calculate intermutation distance
mut_dist_vcfs = lapply(vcfs, function(x) mut_dist(x))

#3 Generate table to check double nucleotides manually
#3A Extract double nucleotides
confirm <- mut_dist_vcfs$`Ercc1(-/D)1 Liver`[which(mut_dist_vcfs$`Ercc1(-/D)1 Liver`[,3] == 1),]
confirm <- rbind(confirm,mut_dist_vcfs$`Ercc1(-/D)1 Small intestine`[which(mut_dist_vcfs$`Ercc1(-/D)1 Small intestine`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`Ercc1(-/D)2 Liver`[which(mut_dist_vcfs$`Ercc1(-/D)2 Liver`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`Ercc1(-/D)2 Small intestine`[which(mut_dist_vcfs$`Ercc1(-/D)2 Small intestine`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`Ercc1(-/D)3 Liver`[which(mut_dist_vcfs$`Ercc1(-/D)3 Liver`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`Ercc1(-/D)3 Small intestine`[which(mut_dist_vcfs$`Ercc1(-/D)3 Small intestine`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`WT1 Liver`[which(mut_dist_vcfs$`WT1 Liver`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`WT2 Liver`[which(mut_dist_vcfs$`WT2 Liver`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`WT2 Small intestine`[which(mut_dist_vcfs$`WT2 Small intestine`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`WT3 Liver`[which(mut_dist_vcfs$`WT3 Liver`[,3] == 1),])
confirm <- rbind(confirm,mut_dist_vcfs$`WT3 Small intestine`[which(mut_dist_vcfs$`WT3 Small intestine`[,3] == 1),])
#3B Put in a table and add sample type
dinucleotides <- data.frame(
  chr = confirm$chromosome,
  pos = confirm$location-1,
  sample <- c(rep("Ercc1(-/D)1 Liver",nrow(mut_dist_vcfs$`Ercc1(-/D)1 Liver`[which(mut_dist_vcfs$`Ercc1(-/D)1 Liver`[,3] == 1),])),
              rep("Ercc1(-/D)1 Small intestine",nrow(mut_dist_vcfs$`Ercc1(-/D)1 Small intestine`[which(mut_dist_vcfs$`Ercc1(-/D)1 Small intestine`[,3] == 1),])),
              rep("Ercc1(-/D)2 Liver",nrow(mut_dist_vcfs$`Ercc1(-/D)2 Liver`[which(mut_dist_vcfs$`Ercc1(-/D)2 Liver`[,3] == 1),])),
              rep("Ercc1(-/D)2 Small intestine",nrow(mut_dist_vcfs$`Ercc1(-/D)2 Small intestine`[which(mut_dist_vcfs$`Ercc1(-/D)2 Small intestine`[,3] == 1),])),
              rep("Ercc1(-/D)3 Liver",nrow(mut_dist_vcfs$`Ercc1(-/D)3 Liver`[which(mut_dist_vcfs$`Ercc1(-/D)3 Liver`[,3] == 1),])),
              rep("Ercc1(-/D)3 Small intestine",nrow(mut_dist_vcfs$`Ercc1(-/D)3 Small intestine`[which(mut_dist_vcfs$`Ercc1(-/D)3 Small intestine`[,3] == 1),])),
              rep("WT1 Liver",nrow(mut_dist_vcfs$`WT1 Liver`[which(mut_dist_vcfs$`WT1 Liver`[,3] == 1),])),
              rep("WT2 Liver",nrow(mut_dist_vcfs$`WT2 Liver`[which(mut_dist_vcfs$`WT2 Liver`[,3] == 1),])),
              rep("WT2 Small intestine",nrow(mut_dist_vcfs$`WT2 Small intestine`[which(mut_dist_vcfs$`WT2 Small intestine`[,3] == 1),])),
              rep("WT3 Liver",nrow(mut_dist_vcfs$`WT3 Liver`[which(mut_dist_vcfs$`WT3 Liver`[,3] == 1),])),
              rep("WT3 Small intestine",nrow(mut_dist_vcfs$`WT3 Small intestine`[which(mut_dist_vcfs$`WT3 Small intestine`[,3] == 1),]))
              )
)
colnames(dinucleotides) <- c("chr","pos","sample")
remove(confirm)
#3C Write table
#write.table(dinucleotides,paste(outdir,"Dinucleotides.txt",sep = ""),sep="\t", row.names = FALSE, col.names = TRUE)
#3D These mutations were checked manually in IGV




# ---- GET DATA: DOUBLE COUNTS ----

#1 Import table: IGV-verfied double base substitutions
checkeddouble <- read.delim2(paste(indir,"SNV/dinucleotides_igv.txt",sep=""))
checkeddouble <- checkeddouble[which(checkeddouble$eyeball == "TP"),]

#2 Import table: Callable genome size
doublecorrectedmutnumber <- read.delim2(paste(indir,"callable_genome/Callable_ERCC1.txt", sep = ""))
doublecorrectedmutnumber$type <- c(rep(c("Ercc1-/D Liver","Ercc1-/D Small intestine"),3),rep(c("WT Liver","WT Small intestine"),3))
doublecorrectedmutnumber$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(doublecorrectedmutnumber$Surveyed.percentage)))

#3 Count double base substitutions & add to callable genome table
doublecorrectedmutnumber$uncorr.double <- c(as.numeric(length(which(checkeddouble$mousetype == "Ercc1(-/D)1 Liver"))),
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

#4 Calculate the corrected number of double base substitutions (extrapolate to autosomal genome)
doublecorrectedmutnumber$corr.double <- (doublecorrectedmutnumber$uncorr.double/doublecorrectedmutnumber$Surveyed.percentage)*100

#5 Add uncorrected number of SNVs to table, to calculate percentage of SNVs that are next to eachother
doublecorrectedmutnumber$uncorr.snv.di <- NA
for (i in vcf_files) {
  df.temp <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
  snv_number <- nrow(df.temp)
  
  # Add to callable genome table
  doublecorrectedmutnumber[which(as.factor(toupper(substr(doublecorrectedmutnumber$Sample,1,5))) == toupper(substr(sample,6,10))),]$uncorr.snv.di <- snv_number
}
remove(df.temp,i,sample,snv_number)
# Subtract Dinucleotide changes (these are in the SNV vcfs)
doublecorrectedmutnumber$uncorr.snv <- doublecorrectedmutnumber$uncorr.snv-(doublecorrectedmutnumber$uncorr.double *2)
# Correct for callable genome (extrapolate SNVs to autosomal genome)
doublecorrectedmutnumber$corr.snv <- (doublecorrectedmutnumber$uncorr.snv/doublecorrectedmutnumber$Surveyed.percentage)*100

#6 Correct number of double base substitutions for total number of SNVs
doublecorrectedmutnumber$corr.double.fraction <- doublecorrectedmutnumber$corr.double/doublecorrectedmutnumber$corr.snv

#7 Calculate mutation rate per week
doublecorrectedmutnumber$corr.double.perweek <- doublecorrectedmutnumber$corr.double/16
#write.table(doublecorrectedmutnumber,file = paste(outdir,"doublesnvnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)

#8 Generate data frame with all 12 samples (Ercc1/WT & Liver/SI)
doublemutnumber <- data.frame(
  type = factor(c("WT_Liver","ERCC1D-_Liver","WT_SmallIntestine","ERCC1D-_SmallIntestine"),
                levels=c("WT_Liver","ERCC1D-_Liver","WT_SmallIntestine", "ERCC1D-_SmallIntestine")),
  mouse = c("WT1 Liver","Ercc1(-/D)1 Liver",
            "WT1 Small intestine","Ercc1(-/D)1 Small intestine",
            "WT2 Liver","Ercc1(-/D)2 Liver",
            "WT2 Small intestine","Ercc1(-/D)2 Small intestine",
            "WT3 Liver","Ercc1(-/D)3 Liver",
            "WT3 Small intestine","Ercc1(-/D)3 Small intestine"), 
  mutations.perweek = doublecorrectedmutnumber[c(7,1,8,2,9,3,10,4,11,5,12,6),]$corr.double.perweek,
  mutations.fractionSNV = doublecorrectedmutnumber[c(7,1,8,2,9,3,10,4,11,5,12,6),]$corr.double.fraction
)




# ---- STATISTICAL ANALYSIS ----

#1 Welch Two Sample t-test
pliver.double <- t.test(doublemutnumber[grep("Liver",doublemutnumber$type),"mutations.perweek"] ~ doublemutnumber[grep("Liver",doublemutnumber$type),"type"])$p.value
psi.double <- t.test(doublemutnumber[grep("Small",doublemutnumber$type),"mutations.perweek"] ~ doublemutnumber[grep("Small",doublemutnumber$type),"type"])$p.value
pliver.doublefraction <- t.test(doublemutnumber[grep("Liver",doublemutnumber$type),"mutations.fractionSNV"] ~ doublemutnumber[grep("Liver",doublemutnumber$type),"type"])$p.value
psi.doublefraction <- t.test(doublemutnumber[grep("Small",doublemutnumber$type),"mutations.fractionSNV"] ~ doublemutnumber[grep("Small",doublemutnumber$type),"type"])$p.value

#2 Adjust for multiple testing
p.adjust(c(pliver.double,psi.double),method = "fdr") #0.002641489 0.999783221
p.adjust(c(pliver.doublefraction,psi.doublefraction), method = "fdr") #0.02862091 0.85404197

#3 Poisson test
pois.liver <- poisson.test(x = c(round(2*sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "Ercc1-/D Liver"),"corr.double"])), 
                   round(sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "Ercc1-/D Liver"),"corr.snv"]))),
             T = c(round(2*sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "WT Liver"),"corr.double"])), 
                   round(sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "WT Liver"),"corr.snv"])))
)
p.pois.liver<- pois.liver$p.value

pois.si <- poisson.test(x = c(round(2*sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "Ercc1-/D Small intestine"),"corr.double"])), 
                   round(sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "Ercc1-/D Small intestine"),"corr.snv"]))),
            T = c(round(2*sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "WT Small intestine" & doublecorrectedmutnumber$Sample != "WT1SI"),"corr.double"])), 
                   round(sum(doublecorrectedmutnumber[which(doublecorrectedmutnumber$type == "WT Small intestine"& doublecorrectedmutnumber$Sample != "WT1SI"),"corr.snv"])))
)
p.pois.si <- pois.si$p.value

#4 Adjust for multiple testing
p.adjust(c(p.pois.liver,p.pois.si), method = "fdr") # 4.974417e-18 1.000000e+00




# ---- PLOT ----

#1 Per ASC type: get statistics (for plotting)
mutdouble.perweek = summarySE(doublemutnumber[-3,], measurevar = "mutations.perweek", groupvars = c("type"))
mutdouble.fractionSNV = summarySE(doublemutnumber[-3,], measurevar = "mutations.fractionSNV", groupvars = c("type"))

#2 Calculate max and min values (for plotting)
mutdouble.perweek <- remove.negative.sd(mutdouble.perweek)
mutdouble.fractionSNV <- remove.negative.sd(mutdouble.fractionSNV)

#3A Plot double point mutation counts per ASC type
doublenumberplot <- ggplot(data = mutdouble.perweek, aes(x=type,y=mutations.perweek,fill=type)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin=neg, ymax=pos), width=.07, position=position_dodge(.9)) +
  xlab(" Ercc1") + 
  ylab("Double base substitutions\nper autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=0.61) +
  ggtitle("") +
  guides(fill=FALSE) +  
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  scale_fill_manual(values=c("lightcoral","skyblue3",
                             "lightcoral","skyblue3")) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())

#3B Save plot
#ggsave(paste(outdir,"DOUBLEnumber.pdf",sep=""), doublenumberplot , width = 8, height = 8)

#4A Plot double point mutation counts per ASC
doublenumberpersampleplot <- ggplot(data = doublemutnumber, aes(x=type,y=mutations.perweek,fill=mouse)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab(" Ercc1") + 
  ylab("Double base substitutions\nper autosomal genome per week") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=0.61) +
  guides(fill=FALSE) +  #geen legenda
  ggtitle("") +
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  scale_fill_manual(values=c("skyblue3","skyblue3","lightskyblue2","lightskyblue2","steelblue3","steelblue3",
                             "lightpink","lightpink","lightcoral","lightcoral","indianred3","indianred3")) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())

#4B Save plot
#ggsave(paste(outdir,"DOUBLEnumber_permouse.pdf",sep=""), doublenumberpersampleplot, width = 8, height = 8)

#5A Plot fraction double base substitutions of total base substitutions per ASC
doublenumberplot.fraction <- ggplot(data = doublemutnumber, aes(x=type,y=mutations.fractionSNV,fill=mouse)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab(" Ercc1") + 
  ylab("Fraction\n(double base substitutions/base substitutions)\n") +
  scale_y_continuous(expand = c(0, 0)) +
 # expand_limits(y=0.002) +
  ggtitle("") +
  guides(fill=FALSE) +  
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  scale_fill_manual(values=c("skyblue3","skyblue3","lightskyblue2","lightskyblue2","steelblue3","steelblue3",
                             "lightpink","lightpink","lightcoral","lightcoral","indianred3","indianred3")) + 
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())

#5B Save plot
#ggsave(paste(outdir,"DOUBLEnumber_permouse_fraction.pdf",sep=""), doublenumberplot.fraction, width = 8, height = 8)
