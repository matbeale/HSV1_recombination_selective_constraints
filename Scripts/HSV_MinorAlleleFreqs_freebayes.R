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
if(!require(gridExtra)){
  install.packages("gridExtra",repos="http://www.stats.bris.ac.uk/R/")
  library("gridExtra")
}
if(!require(grid)){
  install.packages("grid",repos="http://www.stats.bris.ac.uk/R/")
  library("grid")
}
#if(!require(ggExtra)){
#  install.packages("ggExtra",repos="http://www.stats.bris.ac.uk/R/")
#  library("ggExtra")
#}
#if(!require(vcfR)){
 # install.packages("vcfR",repos="http://www.stats.bris.ac.uk/R/")
  #library(vcfR)
#}
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
#biocLite("GenomicFeatures")
#biocLite("BSgenome")
library("VariantAnnotation")
library("GenomicFeatures")
library("BSgenome")


args <- commandArgs(trailingOnly = TRUE) # capture command line arguments
# Put command line arguments into variables
vcf.file <- args[1]

#vcf.file <-  "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/bams/HSV1-CSF1.freebayes.5-reads.ref.vcf"
seqpath <- sub("HSV1-.+","",vcf.file)
seqname <- sub("^.+HSV1-","HSV1-",vcf.file,perl=T)
seqname <- sub("\\..+$","",seqname,perl=T)


# create transcript database from gff - may want to tweak this later
#hsv1.db <- makeTxDbFromGFF(file="C:/Bioinformatics/HSV1_genomes/annotation/NC_001806.2.gff",format="gff3",dataSource="NCBI")
#saveDb(hsv1.db,file="C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/hsv1.db.sqlite")
#txdb <- loadDb("C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/hsv1.db.sqlite")

#refgenome.file <- "C:/Bioinformatics/HSV1_genomes/annotation/NC_001806.2.fa"
#refgenome <- FaFile(refgenome.file)


#intersect(seqlevels(vcf), seqlevels(txdb)) # double check levels (naming) of vcf and gff-db
#codingloc <- locateVariants(vcf, txdb, CodingVariants())
#alternate <- as.character(unlist(alt(vcf)))[values(codingloc)[["QUERYID"]]]
#reference <- as.character(ref(vcf))[values(codingloc)[["QUERYID"]]]
#annotable <- cbind(as.data.frame(codingloc,row.names = NULL), alternate, reference)#[-c(5,7,8,9)]

#coding <- predictCoding(vcf, txdb, seqSource=refgenome)
#####
#coding.df <- as.data.frame( coding, "data.frame", row.names = NULL)
#coding.df <- coding.df[,c(-6,-8)]
#coding.df$CDSID <- gsub("c\\(","" , coding.df$CDSID)
#coding.df$CDSID <- gsub("\\)","" , coding.df$CDSID)
#coding.df$CDSID <- gsub("\\, ","," , coding.df$CDSID)

#myvcf <- data.frame(fixed(vcf))

#myvcf <- data.frame(info(vcf))
#myvcf$position <- as.numeric(gsub("\\_.+$","",gsub("^.+\\:","", rownames(myvcf),perl=T), perl=T))


#working.folder <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/"
#seq.names <- gsub(".varscan.subsampled-500000.txt","",list.files(path=varscan.folder, pattern="subsampled-500000"))

#varscan2file <- paste(working.folder,currentsample,".varscan.5-reads.txt", sep="")

#snps.out = NULL

#vcf.file <- "C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/MinorVars/HSV1-CSF1.freebayes-par.5-reads.vcf"


#########################################################  
## Extract vcf into dataframe
##########################

# with vcfR
#myvcf <- vcfR2tidy(read.vcfR(vcf.file),single_frame=T)
#myvcf.fix <- myvcf$dat # put all data (except meta) into df
#myvcf.pruned <- myvcf.fix[,c(1,2,4,5,8,9,12,14,15,22,23,24,25)]
#colnames(myvcf.pruned)[c(6:9)] <- c("Depth","Alleles","RefCount","AltCount")

# Do this using VariantAnnotation package
#vcf.header <- scanVcfHeader(vcf.file)
vcf <- readVcf(vcf.file,"NC_001806.2_HSV1_s17")
myvcf <- data.frame(info(vcf))
myvcf$position <- as.numeric(gsub("\\_.+$","",gsub("^.+\\:","", rownames(myvcf),perl=T), perl=T))



myvcf.pruned <- myvcf[,c("position","NS","DP","AN","RO","AO","SRF","SRR","SAF","SAR")]
colnames(myvcf.pruned)[c(3:6)] <- c("Depth","Alleles","RefCount","AltCount")
################################

# Seems like the subsampling/depth reduction produces an increase in complex variant calls in FreeBayes - need to be careful
# Need to filter out multi-frequency variant positions (complex variants with weird frequencies - messes up calculations)
myvcf.pruned$keep <- sapply(1:nrow(myvcf.pruned), function(x) ifelse((grepl(",",myvcf.pruned$AltCount[x]) ),"drop","keep"))
myvcf.pruned$keep <- sapply(1:nrow(myvcf.pruned), function(x) ifelse((grepl(",",myvcf.pruned$RefCount[x]) ),"drop",myvcf.pruned$keep[x]))
myvcf.pruned <- subset(myvcf.pruned,keep=="keep", select=c(1:ncol(myvcf.pruned)-1))
# Now resume analysis
myvcf.pruned$Depth2 <- sapply(1:nrow(myvcf.pruned), function(x) sum(as.numeric(myvcf.pruned$RefCount[x]), as.numeric(myvcf.pruned$AltCount[x]))) # depth provided in VCF is missing some reads, so recalculate this
myvcf.pruned$refper <- round(((as.numeric(myvcf.pruned$RefCount)/as.numeric(myvcf.pruned$Depth2))*100),1)
myvcf.pruned$altper <- round(((as.numeric(myvcf.pruned$AltCount)/as.numeric(myvcf.pruned$Depth2))*100),1)

# Produce filter to remove "major variant sites"
myvcf.pruned$vartype <- sapply(1:nrow(myvcf.pruned), function(x) ifelse(myvcf.pruned$refper[x] >= 98 || myvcf.pruned$altper[x] >=98, "major","minor" ))
myvcf.minors <- subset(myvcf.pruned, vartype=="minor")
# Modify filter to remove minor vars with low coverage
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$RefCount[x]) <5, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$AltCount[x]) <5, "drop",myvcf.minors$vartype[x]))
myvcf.minors <- subset(myvcf.minors, vartype=="minor")
myvcf.minors$majorfreq <- sapply(1:nrow(myvcf.minors), function (x) max(c(as.numeric(myvcf.minors$refper[x]) ,as.numeric(myvcf.minors$altper[x]))))
myvcf.minors$minorfreq <- sapply(1:nrow(myvcf.minors), function (x) min(c(as.numeric(myvcf.minors$refper[x]) ,as.numeric(myvcf.minors$altper[x]))))
  
# Add in filter to remove minor variants not covered by reads in both strands
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SRF[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SRR[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SAF[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SAR[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors <- subset(myvcf.minors, vartype=="minor")


if(nrow(myvcf.minors)==0) {
  print("no minor variants in subsampled freebayes data")
  write.table(file=paste(seqpath,seqname,".MinorVars.freebayes.txt",sep=""),myvcf.minors, col.names=T, sep="\t", quote=T)
}



if(nrow(myvcf.minors)!=0) {
# plot data
genomelength <- 152222
windowsize <- 100 #specify size of window
binsize <- 2

# Put all minor variants into single value (i.e. vs major)
wgs.allele.freqs <- data.frame(c(1:152222))
colnames(wgs.allele.freqs) <- "Position"
wgs.allele.freqs <- merge(wgs.allele.freqs,myvcf.minors, by.x="Position", by.y="position", all.x=T)
wgs.allele.freqs$minorfreq <- as.numeric(wgs.allele.freqs$minorfreq)
nwindows <- nrow(wgs.allele.freqs)/windowsize #determine number of windows in subset
wgs.allele.freqs$window <- (trunc(wgs.allele.freqs$Position / windowsize,0))*windowsize # create binning variables and re-label to nearest bin 
wgs.allele.freqs2 <- wgs.allele.freqs # create duplicate df for downstream work before subsetting
wgs.allele.freqs <- subset(wgs.allele.freqs, vartype=="minor")
coverage.window <- aggregate(wgs.allele.freqs$minorfreq, list(wgs.allele.freqs$window), max) # aggregate depth based on bin and calculate mean
names(coverage.window) <- c("Position", "Max_Vars")
# Plot minor variants along genome
p2 <- ggplot(coverage.window, aes(Position, Max_Vars))
p2 <- p2 + geom_bar(stat="identity", colour="steelblue4") +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  labs(x=paste("Genome Position (bp; windowsize=",windowsize,")",sep=""), y="Minor Variant Frequency (%)") +
  scale_x_continuous(labels = scales::comma,expand = c(0, 0), breaks=seq(0, genomelength, round((genomelength / 15),1-nchar(ceiling(genomelength[1] / 15)))),limits = c(0, genomelength)) +
  scale_y_continuous(expand=c(0, 0),limits = c(0, 50+binsize)) + # add binsize from frequency plot so that y axes line up
  annotate("text", x = 12000, size=5, y = 48, label = seqname)
  #ggtitle(seqname) 
p2
# Now do frequency plot
#binsize <- 2
wgs.allele.freqs$freqbin <- ((ceiling(wgs.allele.freqs$minorfreq / binsize))*binsize)-(binsize/2) # create binning variables and re-label to nearest bin 
p3 <- ggplot(wgs.allele.freqs, aes(freqbin))
p3 <- p3 + geom_bar(stat="count", colour="steelblue4", fill="steelblue3") +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  labs(x="", y="Occurrences in Genome (count)") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous( expand=c(0,0),breaks=seq(0,50,binsize),limits = c(0, 50+binsize)) + # #need to specify binsize in limits or 48-50 is excluded
  theme(axis.text.x = element_text(hjust = 1, vjust=0.5)) +
  coord_flip() #+
  #theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())
p3
gA <- ggplotGrob(p2 + theme(legend.justification=c(0,0), legend.position=c(0,0) ,legend.text = element_text(size = 10 ),legend.key.size = unit(0.5,"cm"),legend.direction="horizontal" ))
gB <- ggplotGrob(p3 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank() ))
gB$heights <- gA$heights
grid.newpage()

Cairo(file= paste(seqpath,seqname,".MinorVars.freebayes.png",sep=""), width = 1200, height = 300,type="png",dpi=600, pointsize = 12*600/72,units = "pt")
grid.arrange(gA, gB, nrow = 1, widths=c(5,1))
dev.off()

write.table(file=paste(seqpath,seqname,".MinorVars.freebayes.txt",sep=""),myvcf.minors, col.names=T, sep="\t", quote=T)

}


