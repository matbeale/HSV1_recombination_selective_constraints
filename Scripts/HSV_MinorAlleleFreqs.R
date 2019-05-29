#!/usr/bin/env Rscript

#Usage : HSV_MinorAlleleFreqs.R \
#  <Consensus.fasta>   # 1
#  <varscan2file.txt>  # 2

args <- commandArgs(trailingOnly = TRUE) # capture command line arguments
# Put command line arguments into variables
consensusfile <- args[1]
varscan2file <- args[2]

#consensusfile <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/denovo/Sofias/Consensus/HSV-Matt-009.consensus.fa"
#varscan2file  <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/denovo/Sofias/Minor_variants/uncorrected/HSV-Matt-009.minorvars"


# dependencies
if(!require(Cairo)){
  install.packages("Cairo",repos="http://www.stats.bris.ac.uk/R/")
  library("Cairo")
}
if(!require(ggplot2)){
  install.packages("ggplot2",repos="http://www.stats.bris.ac.uk/R/")
  library(ggplot2)
}
if(!require(seqinr)){
  install.packages("seqinr",repos="http://www.stats.bris.ac.uk/R/")
  library(seqinr)
}
if(!require(reshape2)){
  install.packages("reshape2",repos="http://www.stats.bris.ac.uk/R/")
  library(reshape2)
}
if(!require(gridExtra)){
  install.packages("gridExtra",repos="http://www.stats.bris.ac.uk/R/")
  library(gridExtra)
}
#if(!require(grid)){
#  install.packages("grid",repos="http://www.stats.bris.ac.uk/R/")
#  library("grid")
#}

# Capture file names & paths
seqname <- gsub("^.*\\/","" ,consensusfile, perl=T)
seqname <- gsub("\\..*$", "", seqname, perl=T)
seqpath <- gsub("\\.varscan.*$", "", varscan2file, perl=T)

# Pull in data
wgsconsensus <-(read.fasta(file=consensusfile,as.string=F, seqtype=c("DNA"),seqonly=F,set.attributes=F))
names(wgsconsensus) <-gsub("\\..*$","",names(wgsconsensus), ignore.case=TRUE, perl=T) # strip down sequence names
wgscols <- data.frame(unlist(wgsconsensus,use.names=F,recursive=F),stringsAsFactors=F)  # take all seqs from list into dataframe columnsusing lapply style list loop
colnames(wgscols) <- names(wgsconsensus) # rename columns in new dataframe
wgscols$position <- c(1:nrow(wgscols))

# Process Varscan data, including dealing with tri-allelic sites
varscan2 <- read.table(varscan2file,header=T, stringsAsFactors=F)
varscan2$VarFreq <- as.numeric(gsub("%","",varscan2$VarFreq))

if(is.data.frame(varscan2) && nrow(varscan2)==0) {
  varscan2m <- rbind(c(".",1,rep(".",17)))
  colnames(varscan2m) <- colnames(varscan2)
  varscan2 <- data.frame(varscan2m, stringsAsFactors=F)}
#if(is.data.frame(varscan2) && nrow(varscan2)!=0) {
varscan2 <- varscan2[,c(2:3,19,7)]
varscan2$Ref <- tolower(varscan2$Ref)
varscan2$VarAllele <- tolower(varscan2$VarAllele)
# create new merged columns (VarAllelleM1,M2, and VarFreq) for duplicate positions (20%)
# Create first minor allele (Var2)
combined <- varscan2$VarAllele[1]
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarAllele[x-1]), paste(varscan2$VarAllele[x]))),after = 1)}
varscan2$VarAlleleM2 <- combined
# Create second minor allele (Var3)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[2], paste(varscan2$VarAllele[1]), paste("-")) 
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarAllele[x]), paste("-"))),after = 1)}
varscan2$VarAlleleM3 <- combined
# Create third minor allele (Var4)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[3], paste(varscan2$Position[1]), c("-")) 
combined <- if(varscan2$Position[1] != varscan2$Position[3]) {append(combined,"-",after=1)}
if(nrow(varscan2)>2) {combined <- append(combined, sapply(3:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-2], paste(varscan2$VarAllele[x]), paste("-"))),after = 1)}
varscan2$VarAlleleM4 <- combined

# Now do the same for Allele Frequencies
combined <- varscan2$VarFreq[1]
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarFreq[x-1]), paste(varscan2$VarFreq[x]))),after = 1)}
varscan2$VarFreqM2 <- combined
# Create second minor allele (Var3)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[2], paste(varscan2$VarFreq[1]), paste("0")) 
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarFreq[x]), paste("0"))),after = 1)}
varscan2$VarFreqM3 <- combined
# Create third minor allele (Var4)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[3], paste(varscan2$Position[1]), c("0")) 

#combined <- if(varscan2$Position[1] != varscan2$Position[3]) {append(combined,"0",after=1)}
if(nrow(varscan2)>2) {combined <- if(varscan2$Position[1] != varscan2$Position[3]) {append(combined,"0",after=1)}}
if(nrow(varscan2)>2) {combined <- append(combined, sapply(3:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-2], paste(varscan2$VarFreq[x]), paste("0"))),after = 1)}
varscan2$VarFreqM4 <- combined

# work out "major variant frequency"
combined <- sapply(1:nrow(varscan2), function(x) 100-((as.numeric(varscan2$VarFreqM2[x]) + as.numeric(varscan2$VarFreqM3[x]) + as.numeric(varscan2$VarFreqM4[x]))))
varscan2$VarFreqM1 <- combined

# Label "extra" (non-merged) duplicate row to drop (20%)
if(nrow(varscan2)>1) {combined <- sapply(1:(nrow(varscan2)-1), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x+1], "drop", "keep"))
combined <- append(combined,"keep", length(combined))
varscan2$keep <- combined}
# subset based on drop (20%)
if(nrow(varscan2)>1) {varscan2 <- subset (varscan2, varscan2$keep =="keep")}

# Combine Minor Variants with WGS Sequence info
wgsminor <- merge(wgscols, varscan2, by.x="position", by.y="Position", all=T)
wgsminor[is.na(wgsminor)] <- "."
wgsminor$VarFreqM1 <- ifelse(wgsminor$VarFreqM1 ==".", 100,wgsminor$VarFreqM1)
wgsminor$VarFreqM2 <- ifelse(wgsminor$VarFreqM2 ==".", 0,wgsminor$VarFreqM2)
wgsminor$VarFreqM3 <- ifelse(wgsminor$VarFreqM3 ==".", 0,wgsminor$VarFreqM3)
wgsminor$VarFreqM4 <- ifelse(wgsminor$VarFreqM4 ==".", 0,wgsminor$VarFreqM4)

# Create new dataframe with only frequencies
#wgsallelefreqs <- data.frame(t(sapply(1:nrow(wgsminor), function(x) sort(wgsminor[x,c(9:12)], decreasing=T))), stringsAsFactors = F)
wgs.allele.freqs <- data.frame(t(sapply(1:nrow(wgsminor), function(x) sort(as.numeric(wgsminor[x,c(9:12)]), decreasing=T))), stringsAsFactors = F)
names(wgs.allele.freqs) <- c("VarFreqM1","VarFreqM2","VarFreqM3","VarFreqM4")
# Put position names back in
wgs.allele.freqs$position <- c(1:nrow(wgs.allele.freqs))
# Strip out major variant frequency
wgs.allele.freqs <- wgs.allele.freqs[c(5,2:4)]

# Reshape (melt) data into longform for stacked plot (multiple variants)
wgs.melted <- melt(wgs.allele.freqs, id.vars="position")
names(wgs.melted) <- c("Position", "Allele", "Frequency")



### Use sliding window

windowsize <- 100 #specify size of window

# Put all minor variants into single value (i.e. vs major)
wgs.allele.freqs$allminor <- wgs.allele.freqs$VarFreqM2 + wgs.allele.freqs$VarFreqM3 + wgs.allele.freqs$VarFreqM4  
wgs.allele.freqs$sub_position <- c(1:nrow(wgs.allele.freqs)) # generate subset position column
nwindows <- nrow(wgs.allele.freqs)/windowsize #determine number of windows in subset
wgs.allele.freqs$window <- (trunc(wgs.allele.freqs$position / windowsize,0))*windowsize # create binning variables and re-label to nearest bin 
coverage.window <- aggregate (wgs.allele.freqs$allminor, list(wgs.allele.freqs$window), max) # aggregate depth based on bin and calculate mean
names(coverage.window) <- c("Position", "Max_Vars")


# Plot minor variants along genome
p2 <- ggplot(coverage.window, aes(Position, Max_Vars))
p2 <- p2 + geom_bar(stat="identity", colour="steelblue4") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=paste("Genome Position (bp; windowsize=",windowsize,")",sep=""), y="Minor Variant Frequency (%)") +
  scale_x_continuous(labels = scales::comma,expand = c(0, 0), breaks=seq(0, nrow(wgs.allele.freqs), round(nrow(wgs.allele.freqs) / 15,1 -nchar(ceiling(nrow(wgs.allele.freqs) / 15))))) +
  scale_y_continuous(expand=c(0, 0),limits = c(0, 50)) +
  ggtitle(seqname)
  
#p2


# Look at frequency of minor variants (clusters of specific %s)
allminoralleles <- subset(wgs.melted, Frequency>=2)
allminoralleles$rounded <- round(allminoralleles$Frequency,0)

# bin data
binsize <- 2
#wgs.allele.freqs$window 
#allminoralleles$window <- (trunc(allminoralleles$Frequency / binsize,0))*binsize # create binning variables and re-label to nearest bin 
#allminoralleles$window <- (round(allminoralleles$Frequency / binsize,0))*binsize # create binning variables and re-label to nearest bin 
allminoralleles$window <- ((ceiling(allminoralleles$Frequency / binsize))*binsize)-(binsize/2) # create binning variables and re-label to nearest bin 


## Plot allele counts
p3 <- ggplot(allminoralleles, aes(window))
p3 <- p3 + geom_bar(stat="count", colour="steelblue4", fill="steelblue3") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="", y="Occurrences in Genome (count)") +
  scale_y_continuous() +
  scale_x_continuous( expand=c(0,0),breaks=seq(0,50,binsize),limits = c(0, 50+binsize)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust=0.5)) +
  coord_flip() #+
  #theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())
#p3


#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +

#grid.arrange(p2, p3, ncol=2, widths=c(4,1))
#gA <- ggplotGrob(p2)
#gB <- ggplotGrob(p3)
#gB$heights <- gA$heights
#grid.arrange(gA, gB, ncol=2, widths=c(4,1))

Cairo(file= paste(seqpath,".MinorVars.png",sep=""), width = 900, height = 450,type="png",dpi=900, pointsize = 12*900/72,units = "pt")
#Cairo(file= paste(seqpath,".MinorVars.png",sep=""), width = 900, height = 450,type="png",dpi=900)
#Cairo(file= paste(varscan2file,".MinorVars.png",sep=""), width = 1200, height = 600,type="png",dpi=600, pointsize = 12*600/72,units = "pt")
#tiff(paste(seqpath,".MinorVars.tiff",sep=""),width = 1200, height = 600)
grid.arrange(p2, p3, ncol=2, widths=c(4,1))
#grid.arrange(gA, gB, ncol=2, widths=c(4,1))
dev.off()
#p4


