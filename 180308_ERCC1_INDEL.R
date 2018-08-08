# @Date: 30 October 2017
# @Modified: 08 March 2018
# @Author: Myrthe Jager
# @Description: Indel analysis ERCC1 delta min mouse organoids
# Abbreviations: MUT = Ercc1 delta min; WT = Wild-type, SI = Small intestine, DEL = deletion, IN = insertion




# ---- GET STARTED ----

#1 Install & load required packages
#install.packages("ggplot2")
library("ggplot2")

#2 Define input and output directory
indir = "~/surfdrive/Shared/ERCC1/Data/" 
outdir = "~/surfdrive/Shared/ERCC1/Results/INDEL/"

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




# ---- GET DATA: INDEL COUNTS ----

#1 Import table: Callable genome
indelcorrectedmutnumber <- read.delim2(paste(indir,"callable_genome/Callable_ERCC1.txt", sep = ""))
indelcorrectedmutnumber$Surveyed.percentage <- as.numeric(gsub("%", "", unlist(indelcorrectedmutnumber$Surveyed.percentage)))

#2 Add empty columns to callable genome table for the counts
indelcorrectedmutnumber$uncorr.INDEL <- NA
indelcorrectedmutnumber$uncorr.IN <- NA
indelcorrectedmutnumber$uncorr.DEL <- NA

#3 List vcf files: filtered indels
indelvcf_files <- list.files(paste(indir,"INDEL/",sep = ""), full.names = T)

#4 Count indels & add to callable genome table
for (i in indelvcf_files) {
  df <- read.table(i)
  sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
  ins <- nrow(df[which(as.character(df$V4) == "A" | as.character(df$V4) == "T" | as.character(df$V4) == "C" | as.character(df$V4) == "G"),]) 
  del <- nrow(df[which(as.character(df$V5) == "A" | as.character(df$V5) == "T" | as.character(df$V5) == "C" | as.character(df$V5) == "G"),]) 
  
  #4B Add to the surveyed file
  indelcorrectedmutnumber[which(as.factor(toupper(substr(indelcorrectedmutnumber$Sample,1,5))) == toupper(substr(sample,6,10))),]$uncorr.INDEL <- ins+del
  indelcorrectedmutnumber[which(as.factor(toupper(substr(indelcorrectedmutnumber$Sample,1,5))) == toupper(substr(sample,6,10))),]$uncorr.IN <- ins
  indelcorrectedmutnumber[which(as.factor(toupper(substr(indelcorrectedmutnumber$Sample,1,5))) == toupper(substr(sample,6,10))),]$uncorr.DEL <- del
  remove(ins, sample, del, df)
}
remove(i)

#5 Calculate the corrected number of INDELS
indelcorrectedmutnumber["corr.INDEL"] <- (indelcorrectedmutnumber$uncorr.INDEL/indelcorrectedmutnumber$Surveyed.percentage)*100
indelcorrectedmutnumber["corr.IN"] <- (indelcorrectedmutnumber$uncorr.IN/indelcorrectedmutnumber$Surveyed.percentage)*100
indelcorrectedmutnumber["corr.DEL"] <- (indelcorrectedmutnumber$uncorr.DEL/indelcorrectedmutnumber$Surveyed.percentage)*100

#6 Calculate mutation rate per week
indelcorrectedmutnumber$corr.INDEL.perweek <- indelcorrectedmutnumber$corr.INDEL/16
indelcorrectedmutnumber$corr.IN.perweek <- indelcorrectedmutnumber$corr.IN/16
indelcorrectedmutnumber$corr.DEL.perweek <- indelcorrectedmutnumber$corr.DEL/16
#write.table(indelcorrectedmutnumber,file = paste(outdir,"indelnumber.txt",sep=""),sep ="\t",col.names = T,row.names = F)

#7 Generate data frame with all 12 samples (Ercc1/WT & Liver/SI)
indelmutnumber <- data.frame(
  type = factor(c("WT Liver","Ercc1-/D Liver","WT Small intestine","Ercc1-/D Small intestine"),
                levels=c("WT Liver","Ercc1-/D Liver","WT Small intestine", "Ercc1-/D Small intestine")),
  mouse = c("WT1 Liver","Ercc1(-/D)1 Liver","WT1 Small intestine","Ercc1(-/D)1 Small intestine",
            "WT2 Liver","Ercc1(-/D)2 Liver","WT2 Small intestine","Ercc1(-/D)2 Small intestine",
            "WT3 Liver","Ercc1(-/D)3 Liver","WT3 Small intestine","Ercc1(-/D)3 Small intestine"), 
  mutationtype = c(rep("indel",12),rep("in",12),rep("del",12)),
  mutations.perweek = c(indelcorrectedmutnumber[c(7,1,8,2,9,3,10,4,11,5,12,6),]$corr.INDEL.perweek
  ,indelcorrectedmutnumber[c(7,1,8,2,9,3,10,4,11,5,12,6),]$corr.IN.perweek
  ,indelcorrectedmutnumber[c(7,1,8,2,9,3,10,4,11,5,12,6),]$corr.DEL.perweek)
)




# ---- STATISTICAL ANALYSIS ----

#1 Welch Two Sample t-test
pliver.indel <- t.test(indelmutnumber[grep("Liver",indelmutnumber$type),][grep("indel",indelmutnumber[grep("Liver",indelmutnumber$type),]$mutationtype),]$mutations.perweek ~ indelmutnumber[grep("Liver",indelmutnumber$type),][grep("indel",indelmutnumber[grep("Liver",indelmutnumber$type),]$mutationtype),]$type)$p.value
psi.indel <- t.test(indelmutnumber[grep("Small",indelmutnumber$type),][grep("indel",indelmutnumber[grep("Small",indelmutnumber$type),]$mutationtype),]$mutations.perweek ~ indelmutnumber[grep("Small",indelmutnumber$type),][grep("indel",indelmutnumber[grep("Small",indelmutnumber$type),]$mutationtype),]$type)$p.value

#2 Adjust for multiple testing
p.adjust(c(pliver.indel,psi.indel), method = "fdr") #0.4252863 0.4252863




# ---- PLOT ----

#1 Per ASC type: get statistics (for plotting)
mutindel.perweek <- summarySE(indelmutnumber[c(-3,-15,-27),], measurevar="mutations.perweek",groupvars=c("type","mutationtype"))
mutindel.perweek.combined <- summarySE(indelmutnumber[c(-3,-15,-27),], measurevar="mutations.perweek",groupvars=c("mutationtype"))

#2 Add counts in and del (for plotting)
mutindel.perweek$mutcount <- mutindel.perweek$mutations.perweek
mutindel.perweek[which(mutindel.perweek$mutationtype == "del"),]$mutcount <- mutindel.perweek[which(mutindel.perweek$mutationtype == "indel"),]$mutcount

#3A Plot indel counts per ASC type
indel_mutnumberplot <- ggplot(data=mutindel.perweek[which(mutindel.perweek$mutationtype != "indel"),], aes(x=type, y=mutations.perweek, fill= mutationtype)) +
  geom_bar(stat="identity", width = 0.7) +
  geom_errorbar(aes(ymin=mutcount-sd, ymax=mutcount+sd),
                width=.07) +
  xlab(" Ercc1") + 
  ylab("Indels\nper autosomal genome per week\n") +
  scale_x_discrete(labels=c("WT\nLiver\n","Ercc1-/D\nLiver\n","WT\nSmall\nintestine","Ercc1-/D\nSmall\nintestine"))+
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=20) +
  guides(fill=FALSE) +
  scale_fill_manual(values=c("#999999","#E69F00")) +
  theme(axis.text=element_text(size = rel(1)), 
        axis.title=element_text(size = rel(2)),
        axis.line.y = element_line(colour = "grey40", size = rel(1)),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey40"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x =element_blank(),
        axis.title.x=element_blank())

#3B Save plot
#ggsave(paste(outdir,"INDELnumber.pdf",sep=""), indel_mutnumberplot, width = 8, height = 8)

#4 Plot indel counts per ASC
indel_mutnumberpersampleplot <- ggplot(data=indelmutnumber[which(indelmutnumber$mutationtype != "indel"),], aes(x=mouse, y=mutations.perweek, fill= mutationtype)) +
  geom_bar(stat="identity", width = 0.7) +
  xlab(" Ercc1") + 
  ylab("Indels\nper autosomal genome per week\n") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("#999999","#E69F00")) +
  expand_limits(y=25) +
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
#ggsave(paste(outdir,"INDELnumber_permouse.pdf",sep=""), indel_mutnumberpersampleplot, width = 8, height = 8)




# ---- ADDITIONAL CALCULATIONS ----
#1 Mean Indel number
mean(indelmutnumber[-c(3,13:36),4])

#2 Standard deviation of indel number
sd(indelmutnumber[-c(3,13:36),4])