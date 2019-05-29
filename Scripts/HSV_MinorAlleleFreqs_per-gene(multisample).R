#!/usr/bin/env Rscript

#Usage : HSV_MinorAlleleFreqs.R \
#  <Consensus.fasta>   # 1
#  <varscan2file.txt>  # 2

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


# get the names of all relevant sequences from the varscan files in a particular directory
varscan.folder <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/"
seq.names <- gsub(".varscan.subsampled-500000.txt","",list.files(path=varscan.folder, pattern="subsampled-500000"))

#opbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/HSV1-CSF1.varscan.5-reads.txt"

# pull in coding information & gene positions
gene.matrix <- read.table("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1_gene-matrix.txt",header =T, stringsAsFactors = F)
gene.matrix.melted <- melt(gene.matrix, id.vars="Position")
gene.matrix.melted$value <- as.factor(gene.matrix.melted$value)
# bring in gene names for text and determine midpoint in gene position for plotting
## Pull in gene positions for mapping hits
gene.pos.file <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified"
gene.pos <- read.table(gene.pos.file, header =F, stringsAsFactors = F, quote="")
colnames(gene.pos) <- c("Start","End","Gene","Strand")
gene.pos$Start <- as.numeric(gene.pos$Start)
gene.pos$End <- as.numeric(gene.pos$End)
gene.pos$midpoint <- as.numeric(gene.pos$Start) + ((as.numeric(gene.pos$End) - as.numeric(gene.pos$Start) )/2)

# create a variable for adding to
counts.out = NULL
wgsminor.compiled = NULL

# now loop through all samples
for (currentsample in seq.names){

  consensusfile <- paste("C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/fastas/",currentsample,".consensus.fa", sep="")
  varscan2file <- paste("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/",currentsample,".varscan.5-reads.txt", sep="")
  
  
# Capture file names & paths
seqname <- gsub("^.*\\/","" ,consensusfile, perl=T)
seqname <- gsub("\\..*$", "", seqname, perl=T)
seqpath <- gsub(paste(seqname,"\\.varscan.*$",sep=""), "", varscan2file, perl=T)

# Pull in data
wgsconsensus <- read.fasta(file=consensusfile ,as.string=F, seqtype=c("DNA"),seqonly=F,set.attributes=F)
names(wgsconsensus) <-gsub("\\..*$","",names(wgsconsensus), ignore.case=TRUE, perl=T) # strip down sequence names
wgscols <- data.frame(unlist(wgsconsensus,use.names=F,recursive=F),stringsAsFactors=F)  # take all seqs from list into dataframe columnsusing lapply style list loop
colnames(wgscols) <- names(wgsconsensus) # rename columns in new dataframe
wgscols$position <- c(1:nrow(wgscols))

# Process Varscan data, including dealing with tri-allelic sites
varscan2 <- read.table(varscan2file,header=T, stringsAsFactors=F)
varscan2$VarFreq <- as.numeric(gsub("%","",varscan2$VarFreq))

if(is.data.frame(varscan2) && nrow(varscan2)==0) { # if the dataframe is empty create a dummy df to prevent downstream errors
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
rm(combined)
# Create second minor allele (Var3)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[2], paste(varscan2$VarAllele[1]), paste("-")) 
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarAllele[x]), paste("-"))),after = 1)}
if (is.na(combined[1])) {combined <- "-"}
varscan2$VarAlleleM3 <- combined
rm(combined)
# Create third minor allele (Var4)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[3], paste(varscan2$Position[1]), paste("-")) 
combined <- if(!is.na(combined)) {if(varscan2$Position[1] != varscan2$Position[3]) {append(combined,"-",after=1)}}
#combined <- if(varscan2$Position[1] != varscan2$Position[3]) {append(combined,"-",after=1)}
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(3:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-2], paste(varscan2$VarAllele[x]), paste("-"))),after = 1)}
if (is.na(combined[1])) {combined <- "-"}
varscan2$VarAlleleM4 <- combined


# Now do the same for Allele Frequencies
combined <- varscan2$VarFreq[1]
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarFreq[x-1]), paste(varscan2$VarFreq[x]))),after = 1)}
varscan2$VarFreqM2 <- combined
# Create second minor allele (Var3)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[2], paste(varscan2$VarFreq[1]), paste("0")) 
if(nrow(varscan2)>=2) {combined <- append(combined, sapply(2:nrow(varscan2), function(x) ifelse(varscan2$Position[x] == varscan2$Position[x-1], paste(varscan2$VarFreq[x]), paste("0"))),after = 1)}
varscan2$VarFreqM3 <- combined
rm(combined)

# Create third minor allele (Var4)
combined <- ifelse(varscan2$Position[1] == varscan2$Position[3], paste(varscan2$Position[1]), c("0")) 
if (is.na(combined[1])) {combined <- "0"}
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


# Subset to remove "minor vars" above 50% 
# - occurs where VarScan disagrees with main pipeline about consensus snp
#varscan2VarFreqM1 <- 
#varscan2 <- subset(varscan2, VarFreqM1 < 2)



# Combine Minor Variants with WGS Sequence info
wgsminor <- merge(wgscols, varscan2, by.x="position", by.y="Position", all=T)
wgsminor[is.na(wgsminor)] <- "."
wgsminor$VarFreqM1 <- ifelse(wgsminor$VarFreqM1 ==".", 100,wgsminor$VarFreqM1)
wgsminor$VarFreqM2 <- ifelse(wgsminor$VarFreqM2 ==".", 0,wgsminor$VarFreqM2)
wgsminor$VarFreqM3 <- ifelse(wgsminor$VarFreqM3 ==".", 0,wgsminor$VarFreqM3)
wgsminor$VarFreqM4 <- ifelse(wgsminor$VarFreqM4 ==".", 0,wgsminor$VarFreqM4)

#head(t(wgsminor$VarFreqM1))

#head(append(wgsminor$VarFreqM1,seqname,after=length(wgsminor$VarFreqM1)))

wgsminor.compiled <- rbind(wgsminor.compiled,t(append(wgsminor$VarFreqM1,seqname,after=length(wgsminor$VarFreqM1))))



##### Look at only the minor variants sites - what genes do they affect? ####

# assign gene positions to snps
minorvar.sites <- subset(wgsminor, keep=="keep")
minorvar.sites$VarFreq <- as.numeric(minorvar.sites$VarFreq)
minorvar.sites <- subset(minorvar.sites, VarFreq < 100)
#minorvar.sites <- merge(minorvar.sites, coding.df, by.x="position", by.y="start")

minorvar.sites$genemappings <- gsub("^, ","",sapply(1:nrow(minorvar.sites), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(minorvar.sites$position[y] >= gene.pos$Start[x] && minorvar.sites$position[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)

#genecounts <- gene.pos$Gene

genecounts <- data.frame(t(sapply(1:nrow(gene.pos), function(x) nrow(subset(minorvar.sites, genemappings==gene.pos$Gene[x])))))
colnames(genecounts) <- gene.pos$Gene
genecounts$sample <- seqname
genecounts <- genecounts[,c(ncol(genecounts),1:ncol(genecounts)-1)]

counts.out<- rbind(counts.out,genecounts)
}

write.table(file=paste(seqpath,"HSV1_All_MinorVars-per-gene.txt"), counts.out,  col.names=T, sep="\t", quote=T)





# see if this can be plotted
counts.melted <- melt(counts.out, id.vars="sample")
counts.melted$cohort <-  sapply(1:nrow(counts.melted), function(x) ifelse(grepl("nCSF", counts.melted$sample[x]),"non-CSF","CSF"))

# remove nCSF7 (known mixed sample/superinfection)
counts.melted <-  subset(counts.melted, sample !="HSV1-nCSF7")

p1 <- ggplot(counts.melted, aes(variable, value, fill=sample, colour=sample))
p1 <- p1 + geom_bar(stat="identity",position="dodge",alpha=1/3) + facet_grid(cohort ~ ., scales="free", space="free", drop=T, margins=F) +
  theme_bw() + labs(x="Gene", y="Minor Variants/Sample") +
  ylim(0,8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),legend.position = "top") 
p1

Cairo(file= paste(seqpath,"HSV1_All_MinorVars-per-gene-bar.png",sep=""), width = 900, height = 400,type="png",dpi=600, units = "pt")
p1
dev.off()

