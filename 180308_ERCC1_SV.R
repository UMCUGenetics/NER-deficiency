# @Date: 30 October 2017
# @Modified: 08 March 2018
# @Author: Myrthe Jager
# @Description: SV analysis ERCC1 delta min mouse organoids
# Abbreviations: MUT = Ercc1 delta min; WT = Wild-type, LIV = Liver, SI = Small intestine, DEL = deletion, DUP = duplication, INV = inversion, TRA = translocation, INS = insertion, SV = structural variation




# ---- GET STARTED ----

#1 Install & load required packages
#install.packages("ggplot2")
library("ggplot2")

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/SV/"

#3 Functions
### Function 1: PER TYPE: STATISTICS ###
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

### Function 2: Calculate whiskers based on sd ###
calculate.range.sd.sv <- function(df) {
  df$pos <- df[,4] + df$sd
  df$neg <- df[,4] - df$sd
  if (length(df[which(df$neg < 0),]$neg) >0) {
    df[which(df$neg < 0),]$neg <- 0
  }
  return(df)
}




# ---- GET DATA: SV COUNTS ----

#1 Import table: IGV-verfied structural variations
structuralvariations <- read.table(paste(indir, "SV/SV_IGV_verified.txt", sep=""), sep = "\t", header = T)

#2 Import table: Callable genome size
svcorrectedmutnumber <- read.delim2(paste(indir,"callable_genome/Callable_ERCC1.txt", sep = ""))
svcorrectedmutnumber$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(svcorrectedmutnumber$Surveyed.percentage)))

#3 Count deletions & add to callable genome table
svcorrectedmutnumber$uncorr.DEL <- c(as.numeric(length(which(structuralvariations$Mouse == "Ercc1-/D1" & structuralvariations$Tissue == "Liver" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "Ercc1-/D1" & structuralvariations$Tissue == "SI" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "Ercc1-/D2" & structuralvariations$Tissue == "Liver" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "Ercc1-/D2" & structuralvariations$Tissue == "SI" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "Ercc1-/D3" & structuralvariations$Tissue == "Liver" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "Ercc1-/D3" & structuralvariations$Tissue == "SI" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "WT1" & structuralvariations$Tissue == "Liver" & structuralvariations$Type == "deletion"))),
                                     NA,
                                     as.numeric(length(which(structuralvariations$Mouse == "WT2" & structuralvariations$Tissue == "Liver" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "WT2" & structuralvariations$Tissue == "SI" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "WT3" & structuralvariations$Tissue == "Liver" & structuralvariations$Type == "deletion"))),
                                     as.numeric(length(which(structuralvariations$Mouse == "WT3" & structuralvariations$Tissue == "SI" & structuralvariations$Type == "deletion"))))

#4 Count other structural variations
sum(structuralvariations$Type != "deletion")
# None

#5 Calculate the corrected number of SVs, extrapolated to the autosomal genome
svcorrectedmutnumber["corr.DEL"] <- (svcorrectedmutnumber$uncorr.DEL/svcorrectedmutnumber$Surveyed.percentage)*100

#6 Calculate mutation rate per week
svcorrectedmutnumber$corr.DEL.perweek <- svcorrectedmutnumber$corr.DEL/16
#write.table(svcorrectedmutnumber,file = paste(outdir,"svnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)

#7 Calculate median length deletions
median(structuralvariations$Size.bp.)

#8 Generate data frame with all 12 samples (Ercc1/WT & Liver/SI)
svmutnumber <- data.frame(
  type = factor(c("WT_Liver","ERCC1D-_Liver","WT_SmallIntestine","ERCC1D-_SmallIntestine"),
                levels=c("WT_Liver","ERCC1D-_Liver","WT_SmallIntestine", "ERCC1D-_SmallIntestine")),
  mouse = c("WT1 Liver","Ercc1(-/D)1 Liver",
            "WT1 Small intestine","Ercc1(-/D)1 Small intestine",
            "WT2 Liver","Ercc1(-/D)2 Liver",
            "WT2 Small intestine","Ercc1(-/D)2 Small intestine",
            "WT3 Liver","Ercc1(-/D)3 Liver",
            "WT3 Small intestine","Ercc1(-/D)3 Small intestine"), 
  mutationtype = c(rep("del",12)),
  mutations.perweek= c(svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "WT1LiverA"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "MUT1Liver"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "WT1SI"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "MUT1SI"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "WT2LiverB"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "MUT2LiverC"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "WT2SI"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "MUT2SI"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "WT3Liver"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "MUT3Liver"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "WT3SI"),"corr.DEL.perweek"],
               svcorrectedmutnumber[which(svcorrectedmutnumber$Sample == "MUT3SI"),"corr.DEL.perweek"]))




# ---- STATISTICAL ANALYSIS ----

#1 Welch Two Sample t-test
pliver.sv <- t.test(svmutnumber[grep("Liver",svmutnumber$type),"mutations.perweek"] ~ svmutnumber[grep("Liver",svmutnumber$type),"type"])$p.value
psi.sv <- t.test(svmutnumber[grep("Small",svmutnumber$type),"mutations.perweek"] ~ svmutnumber[grep("Small",svmutnumber$type),"type"])$p.value

#2 Adjust for multiple testing
p.adjust(c(pliver.sv,psi.sv), method = "fdr") #0.5537142 0.5537142




# ---- PLOT ----

#1 Per ASC type: get statistics (for plotting)
mutsv.perweek <- summarySE(svmutnumber[c(-3),], measurevar="mutations.perweek",groupvars=c("type","mutationtype"))

#2 Calculate max and min values (for plotting)
mutsv.perweek <- calculate.range.sd.sv(mutsv.perweek)

#3A Plot SV counts per ASC type
svnumberplot <- ggplot(data=mutsv.perweek, aes(x=type, y=mutations.perweek, fill= mutationtype)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin=neg, ymax=pos),
                width=.07,
                position=position_dodge(.9)) +
  xlab(" Ercc1") + 
  ylab("Structural variations\nper autosomal genome per week\n") +
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=0.5) +
  scale_fill_manual(values=c("#999999","#E69F00")) +
  guides(fill=FALSE) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())

#3B Save plot
#ggsave(paste(outdir, "SVnumber.pdf", sep=""), svnumberplot, width = 8, height = 8)

#4A Plot SV counts per ASC
svnumberpersampleplot <- ggplot(data=svmutnumber, aes(x=mouse, y=mutations.perweek, fill= mutationtype)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
  xlab(" Ercc1") + 
  ylab("Structural variations\nper autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("#999999","#E69F00")) +
  expand_limits(y=0.5) +
  guides(fill=FALSE) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())

#4B Save plot
#ggsave(paste(outdir, "SVnumber_permouse.pdf", sep=""), svnumberpersampleplot, width = 8, height = 8)