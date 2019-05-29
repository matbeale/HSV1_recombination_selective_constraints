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
  library("reshape2")
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



#args <- commandArgs(trailingOnly = TRUE) # capture command line arguments
# Put command line arguments into variables
#vcf.file <- args[1]

#vcf.file <-  "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/bams/HSV1-CSF1.freebayes.5-reads.ref.vcf"
#vcf.file <-  "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/bams/HSV1-CSF1-subsampled.freebayes.500000.vcf"
#vcf.file <- "C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/MinorVars/All-HSV1.freebayes.vcf"


#seqpath <- sub("HSV1-.+","",vcf.file)
#seqname <- sub("^.+HSV1-","HSV1-",vcf.file,perl=T)
#seqname <- sub("\\..+$","",seqname,perl=T)

seqlist <- c("HSV1-CSF10","HSV1-CSF12","HSV1-CSF13","HSV1-CSF1","HSV1-CSF2","HSV1-CSF3","HSV1-CSF4","HSV1-CSF5","HSV1-nCSF10","HSV1-nCSF11","HSV1-nCSF13","HSV1-nCSF14","HSV1-nCSF15","HSV1-nCSF1","HSV1-nCSF4","HSV1-nCSF5","HSV1-nCSF6","HSV1-nCSF7")

myvar.collating <- data.frame(c(1:152222),stringsAsFactors = F)
colnames(myvar.collating) <- "position"

for (currentseq in seqlist) {
  #seqpath <- sub("HSV1-.+","",vcf.file)
  #seqname <- sub("^.+HSV1-","HSV1-",vcf.file,perl=T)
  #seqname <- sub("\\..+$","",seqname,perl=T)
  
  vcf <- readVcf(paste("C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/",currentseq,".filt.ann.vcf",sep=""),"NC_001806.2_HSV1_s17")
  myvcf <- data.frame(info(vcf))
  myvcf$position <- as.numeric(gsub("\\_.+$","",gsub("^.+\\:","", rownames(myvcf),perl=T), perl=T))
  myvcf <- myvcf[,c("position","NS","DP","AN","RO","AO","SRF","SRR","SAF","SAR","ANN")]
  colnames(myvcf)[c(3:6)] <- c("Depth","Alleles","RefCount","AltCount")
  myvcf$keep <- sapply(1:nrow(myvcf), function(x) ifelse((grepl(",",myvcf$AltCount[x]) ),"drop","keep"))
  myvcf$keep <- sapply(1:nrow(myvcf), function(x) ifelse((grepl(",",myvcf$RefCount[x]) ),"drop",myvcf$keep[x]))
  myvcf <- subset(myvcf,keep=="keep", select=c(1:ncol(myvcf)-1))
  myvcf$Depth2 <- sapply(1:nrow(myvcf), function(x) sum(as.numeric(myvcf$RefCount[x]), as.numeric(myvcf$AltCount[x]))) # depth provided in VCF is missing some reads, so recalculate this
  myvcf$refper <- round(((as.numeric(myvcf$RefCount)/as.numeric(myvcf$Depth2))*100),1)
  myvcf$altper <- round(((as.numeric(myvcf$AltCount)/as.numeric(myvcf$Depth2))*100),1)
  # Produce filter to keep if >50% alt variant
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(myvcf$altper[x] >=50, "alt","ref" ))
  myvcf <- subset(myvcf, vartype=="alt")
  # Modify filter to remove minor vars with low coverage
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$AltCount[x]) <5, "drop",myvcf$vartype[x]))
  myvcf <- subset(myvcf, vartype=="alt")
  # Add in filter to remove variants not covered by at least one read on each strand
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$SAF[x]) <1, "drop",myvcf$vartype[x]))
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$SAR[x]) <1, "drop",myvcf$vartype[x]))
  myvcf <- subset(myvcf, vartype=="alt")
  
  myvar.called <- data.frame(myvcf[,c("position","altper")],stringsAsFactors = F)
  colnames(myvar.called)[ncol(myvar.called)] <- currentseq
  myvar.collating <- merge(myvar.collating,myvar.called, by.x="position",by.y="position",all.x=T)
  myvar.collating[ncol(myvar.collating)] <- sapply(1:nrow(myvar.collating), function(x) ifelse(is.na(myvar.collating[x,ncol(myvar.collating)]),0,1))
}

####
# Write to table and/or read in to avoid above loop (slow)
seqpath <- "C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/"
#write.table(myvar.collating, file=paste(seqpath,"HSV1-All-MajorVar-sites_merged.freebayes.txt", sep=""), quote=T, sep="\t",col.names=T, row.names=F)

myvar.collating <- read.table(paste(seqpath,"HSV1-All-MajorVar-sites_merged.freebayes.txt", sep=""), stringsAsFactors = F, header=T)
colnames(myvar.collating) <- gsub("\\.","\\-",colnames(myvar.collating))

# strip out positions without any variants in dataset
myvar.collating$keep <- sapply(1:nrow(myvar.collating), function(x) sum(myvar.collating[x,2:ncol(myvar.collating)]))
myvar.collating <- subset(myvar.collating, keep>0, select=-keep)
nrow(myvar.collating)

# Look only at minor-var sites occurring in multiple genomes
#myvar.collating2 <- myvar.collating
#myvar.collating2$multiple <- sapply(1:nrow(myvar.collating2), function(x) sum(myvar.collating2[x,2:ncol(myvar.collating2)]))
#myvar.collating2 <- subset(myvar.collating2, multiple>1, select=-multiple)
#nrow(myvar.collating2)




#######################################################################
##############

##Do some plotting of SNP density##


####################################
##
##
if(!require(ggplot2)){
  install.packages("ggplot2",repos="http://www.stats.bris.ac.uk/R/")
  library(ggplot2)
}
if(!require(grid)){
  install.packages("grid",repos="http://www.stats.bris.ac.uk/R/")
  library("grid")
}
if(!require(gridExtra)){
  install.packages("gridExtra",repos="http://www.stats.bris.ac.uk/R/")
  library("gridExtra")
}
###############################
# Bring in positional matrix data 
gene.matrix <- read.table("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1_gene-matrix.txt",header =T, stringsAsFactors = F)
library(reshape2)
gene.matrix.melted <- melt(gene.matrix, id.vars="Position")
gene.matrix.melted$value <- as.factor(gene.matrix.melted$value)
# bring in gene names for text and determine midpoint in gene position for plotting
genetext <- gene.pos[,c("midpoint","Strand","Gene")]
genetext$text <- "genelabel"
genetext <- genetext[,c(1,4,3)]
colnames(genetext) <- c("Position","variable","value")
genetext$Position <- as.numeric(genetext$Position)
gene.matrix.melted <- rbind(gene.matrix.melted,genetext)
gene.matrix.melted$Position <- as.numeric(gene.matrix.melted$Position)
# Plot gene positions
p3 <- ggplot(gene.matrix.melted)
p3 <- p3 + geom_point(size=3, shape=15, data=subset(gene.matrix.melted,variable != "genelabel"), aes(x=Position,y=value, group=value, colour=value)) + labs(x="Genome Position (bp)", y="Strand") +
  scale_x_continuous(limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  theme_bw() +   scale_y_discrete() + ylim("-","+","Gene") + theme(legend.position = "none") 
p3 <- p3 + geom_text(data=subset(gene.matrix.melted,variable == "genelabel"), aes(x=Position, y="Gene", label=value), angle = -90, size=2.5)
p3

# Plot snp density along genome
#myvar.melted <- melt(myvar.collating, id.vars="position")
#myvar.melted <- subset(myvar.melted,value==1)
#myvar.melted$variable <- gsub("HSV1-","",myvar.melted$variable)
#myvar.melted$sampletype <- gsub("CSF.+$","CSF",myvar.melted$variable)
#countbinsize <- 500
#p1 <- ggplot(myvar.melted,aes(position,fill=sampletype)) + facet_grid(variable ~ ., scales="free", space="free", drop=T, margins=F)
#p1 <- p1 + geom_histogram(binwidth=countbinsize, alpha=1/2) + theme_bw() +
#  scale_x_continuous(limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
#  ylim(0,20) +
#  xlab(paste("Genome Position (bp; ",countbinsize,"bp sliding windows)",sep="")) + ylab("SNP site density (v.s. Reference NC_001806.2)") +
#  scale_fill_manual(values = c("#009900","red3"))
#p1




windowsize <- 500 #specify size of window
myvar.collating2 <- myvar.collating
nwindows <- nrow(myvar.collating2)/windowsize #determine number of windows in subset
myvar.collating2$window <- (trunc(myvar.collating2$position / windowsize,0))*windowsize # create binning variables and re-label to nearest bin 
library(plyr)
snps.window <- data.frame(unique(myvar.collating2$window),stringsAsFactors = F)
colnames(snps.window) <- "window"
for (sample in c(2:19)) {
  snps.window <- cbind(snps.window,data.frame(aggregate(myvar.collating2[,sample], list(myvar.collating2$window), sum),stringsasfactors=T)[,2]) # aggregate depth based on bin and calculate mean
  colnames(snps.window)[sample] <- colnames(myvar.collating2)[sample]
}

snps.window$mean <- sapply(1:nrow(snps.window), function(x) mean(as.numeric(snps.window[x,c(2:19)])))

snps.window.melted <- melt(snps.window, id.vars="window")
snps.window.melted$variable <- gsub("HSV1-","",snps.window.melted$variable)
snps.window.melted$sampletype <- gsub("CSF.+$","CSF",snps.window.melted$variable)
snps.window.melted$variable <- factor(snps.window.melted$variable, levels=unique(snps.window.melted$variable))

#snps.window.melted$sampletype <- as.factor(snps.window.melted$sampletype)
#snps.window.melted$sampletype <- factor(c("CSF","nCSF","mean"), levels=c("CSF","nCSF","mean"))


p4 <- ggplot(snps.window.melted,aes(window,value,fill=sampletype)) + facet_grid(variable ~ ., scales="free", space="free", drop=T, margins=F)
p4 <- p4 + geom_bar(stat="identity", alpha=1/2) + theme_bw() +
  scale_x_continuous(limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  ylim(0,20) +
  xlab(paste("Genome Position (bp; ",countbinsize,"bp sliding windows)",sep="")) + ylab("SNP site density (v.s. Reference NC_001806.2)") +
  scale_fill_manual(values = c("#009900","grey32","red3"))
p4

# Arrange plots
gA <- ggplotGrob(p4 + theme(legend.position = "none"))
gB <- ggplotGrob(p3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x= element_blank()))
gB$widths <- gA$widths
grid.newpage()
Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_All-Samples_SNP-density-each.png", width = 1200, height = 1200,type="png",dpi=1200, units = "pt")
grid.arrange(gA, gB, nrow = 2, heights=c(14,1))
dev.off()


#
#
#####################################
# Look at No SNPs gene or gene region
#
myvar.collating3 <- myvar.collating
# bring in gene positions and merge
gene.pos.file <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified"
gene.pos <- read.table(gene.pos.file, header =F, stringsAsFactors = F, quote="")
colnames(gene.pos) <- c("Start","End","Gene","Strand")
gene.pos$Start <- as.numeric(gene.pos$Start)
gene.pos$End <- as.numeric(gene.pos$End)
# Assign each variant to gene
myvar.collating3$genemappings <- gsub("^, ","",sapply(1:nrow(myvar.collating3), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(myvar.collating3$position[y] >= gene.pos$Start[x] && myvar.collating3$position[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)
coding.df <- read.table(file="C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/HSV1_MinorVars-freebayes-AllSNPs.coding-data", header=T, stringsAsFactors=F)#, col.names=T)#, row.names=F)
myvar.collating3 <- merge(myvar.collating3, coding.df, by.x="position", by.y="start", all.x=T)

# Determine HSV1 gene region (UL or US)
myvar.collating3$generegion <- sapply(1:nrow(myvar.collating3),function(x) ifelse(grepl("UL",myvar.collating3$genemappings[x]),"UL",""))
myvar.collating3$generegion <- sapply(1:nrow(myvar.collating3),function(x) ifelse(grepl("US",myvar.collating3$genemappings[x]),"US",myvar.collating3$generegion[x]))
nrow(subset(myvar.collating3, generegion=="UL")) # 2090
nrow(subset(myvar.collating3, generegion=="US")) # 311
 # need to determine what number of sites (variant or otherwise) are included in the US and UL regions - more ifficult with overlapping/inverted genes
generegion.proportions <- data.frame(c(1:152222),stringsAsFactors = F)
colnames(generegion.proportions) <- "position"
generegion.proportions$genemappings <- gsub("^, ","",sapply(1:nrow(generegion.proportions), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(generegion.proportions$position[y] >= gene.pos$Start[x] && generegion.proportions$position[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)
generegion.proportions$generegion <- sapply(1:nrow(generegion.proportions),function(x) ifelse(grepl("UL",generegion.proportions$genemappings[x]),"UL",""))
generegion.proportions$generegion <- sapply(1:nrow(generegion.proportions),function(x) ifelse(grepl("US",generegion.proportions$genemappings[x]),"US",generegion.proportions$generegion[x]))
generegion.proportions$generegion <- sapply(1:nrow(generegion.proportions),function(x) ifelse(grepl("RL",generegion.proportions$genemappings[x]),"RL",generegion.proportions$generegion[x]))
generegion.proportions$generegion <- sapply(1:nrow(generegion.proportions),function(x) ifelse(grepl("RS",generegion.proportions$genemappings[x]),"RS",generegion.proportions$generegion[x]))
nrow(subset(generegion.proportions, generegion=="UL")) # 96,595
nrow(subset(generegion.proportions, generegion=="US")) # 10,550
nrow(subset(generegion.proportions, generegion=="RL")) # 6,150
nrow(subset(generegion.proportions, generegion=="RS")) # 7,794

# Number of variant positions per region
nrow(subset(myvar.collating3, generegion=="UL")) / nrow(subset(generegion.proportions, generegion=="UL"))
# UL = 0.0216 variant sites/bp
nrow(subset(myvar.collating3, generegion=="US")) / nrow(subset(generegion.proportions, generegion=="US"))
# US = 0.0295 variant sites/bp

(nrow(myvar.collating) / (nrow(generegion.proportions)-(6150+7794)))*100
# 2.117473% of mapped genome is variant

((nrow(subset(myvar.collating3, generegion=="UL"))+(nrow(subset(myvar.collating3, generegion=="US")))) /(nrow(subset(generegion.proportions, generegion=="UL")) + nrow(subset(generegion.proportions, generegion=="US")))) *100
# 2.240889 % of UL/US positions are polymorphic

myvar.collating3$multigenome <- sapply(1:nrow(myvar.collating3), function(x) sum(myvar.collating3[x,2:19]))
myvar.multigenome <- subset(myvar.collating3, multigenome >=2)
myvar.singlegenome <- subset(myvar.collating3, multigenome ==1)
nrow(myvar.multigenome) # 1593 SNPs found in multiple genomes
nrow(myvar.singlegenome) # 1446 SNPs only found in one genome
(nrow(subset(myvar.singlegenome, generegion!="")) /(nrow(subset(generegion.proportions, generegion=="UL")) + nrow(subset(generegion.proportions, generegion=="US")))) *100
# 1.077% of SNPs in UL/US regions only occur in a single genome
(nrow(subset(myvar.multigenome, generegion!="")) /(nrow(subset(generegion.proportions, generegion=="UL")) + nrow(subset(generegion.proportions, generegion=="US")))) *100
# 1.16% of SNPs in UL/US regions occur in multiple genomes


#
#
#############################
# Do analysis on my HSV1 sequences in heavily trimmed alignment with Szpara 2014 genomes - all relative to s17 reference

#szpara.aln <- seqinr::read.fasta("C:/Bioinformatics/HSV1_genomes/Szpara2014/HSV1_Corededup+Szpara2014_v1_all-gap-pos-stripped+ref.fas",seqtype="DNA", as.string=T,set.attributes=F)
# need to replace this alignment with v2 (contains all seqs, gap stripped, aligned to ref)
szpara.aln <- seqinr::read.fasta("C:/Bioinformatics/HSV1_genomes/Szpara2014/HSV1_Corededup+Szpara2014_v2_all-gap-pos-stripped+ref.fas",seqtype="DNA", as.string=T,set.attributes=F)

szpara.names <- names(szpara.aln)
szpara.aln <- data.frame(t(seqinr::as.matrix.alignment(as.alignment(nb=47,nam=NULL,seq=szpara.aln,com=NULL))), stringsAsFactors = F)
colnames(szpara.aln) <- szpara.names
szpara.aln$position <- rownames(szpara.aln)

# determine if there is a snp at each position
# need to remove gapped positions
szpara.aln$gapped <- sapply(1:nrow(szpara.aln), function(x) grepl("-",paste(unlist(as.character((szpara.aln[x,]))),collapse="")))
nrow(subset(szpara.aln, gapped=="TRUE", select=-gapped)) # 39053 gapped sites
nrow(subset(szpara.aln, gapped=="FALSE", select=-gapped)) # 113169 ungapped sites retained
szpara.aln <- subset(szpara.aln, gapped=="FALSE", select=-gapped)


# filter for positions that only contain 'n' and one allele at same time
szpara.aln$npos <- sapply(1:nrow(szpara.aln), function(x) grep("n",paste(unlist(as.character((szpara.aln[x,]))),collapse="")))
szpara.aln[szpara.aln=="n"] <- ""
# filter for variant positions
szpara.aln$snpsite <- sapply(1:nrow(szpara.aln), function(x) ifelse(nchar(paste(unique(as.character(szpara.aln[x,1:c(ncol(szpara.aln)-2)])),collapse=""))>=2,"snp","nosnp"))
nrow(subset(szpara.aln, snpsite=="snp")) # 5029 SNP sites
nrow(subset(szpara.aln, snpsite=="nosnp")) # 108140 no-SNP sites

# put variant sites into new variable
szpara.variants <- subset(szpara.aln, snpsite=="snp",select=c(-npos,-snpsite))

# Subdivide into separate variables based on cohort
szpara.CSF <- szpara.variants[,grep("-CSF",szpara.names)]
szpara.nCSF <- szpara.variants[,grep("-nCSF",szpara.names)]
szpara.other <- szpara.variants[,grep("CSF",szpara.names, invert=T)]

# Filter each subgroup for variant positions
szpara.CSF$snpsite <- sapply(1:nrow(szpara.CSF), function(x) ifelse(nchar(paste(unique(as.character(szpara.CSF[x,1:c(ncol(szpara.CSF)-2)])),collapse=""))>=2,"snp","nosnp"))
szpara.nCSF$snpsite <- sapply(1:nrow(szpara.nCSF), function(x) ifelse(nchar(paste(unique(as.character(szpara.nCSF[x,1:c(ncol(szpara.nCSF)-2)])),collapse=""))>=2,"snp","nosnp"))
szpara.other$snpsite <- sapply(1:nrow(szpara.other), function(x) ifelse(nchar(paste(unique(as.character(szpara.other[x,1:c(ncol(szpara.other)-2)])),collapse=""))>=2,"snp","nosnp"))


# subset to variant list per group
szpara.CSF.sub <- subset(szpara.CSF, snpsite=="snp")
szpara.nCSF.sub <- subset(szpara.nCSF, snpsite=="snp")
szpara.other.sub <- subset(szpara.other, snpsite=="snp")

szpara.snpbygroup.sub <- list(rownames(szpara.CSF.sub), rownames(szpara.nCSF.sub),rownames(szpara.other.sub))

szpara.list <- (lapply(szpara.snpbygroup.sub, paste,  collapse=" "))#, file="C:/Bioinformatics/HSV1_genomes/Szpara2014/vennlist1.txt")
szpara.list.df <- as.data.frame(do.call(rbind,szpara.list))

write.table(szpara.list.df, file="C:/Bioinformatics/HSV1_genomes/Szpara2014/vennlist1.txt")

#library("VennDiagram")
library("venneuler")

# Filter each subgroup for variant positions
szpara.snpbygroup <- cbind(szpara.CSF$snpsite,szpara.nCSF$snpsite)
szpara.snpbygroup <- cbind(szpara.snpbygroup, szpara.other$snpsite)
szpara.snpbygroup <- data.frame(szpara.snpbygroup)
colnames(szpara.snpbygroup) <- c("CSF","nCSF","other")
rownames(szpara.snpbygroup) <- rownames(szpara.variants)
szpara.snpbygroup$position <- rownames(szpara.snpbygroup)

szpara.snpbygroup.melted <- melt(szpara.snpbygroup, id.vars="position")
szpara.snpbygroup.melted <- subset(szpara.snpbygroup.melted, value=="snp", select=-value)

szpara.vd <- venneuler(szpara.snpbygroup.melted)
szpara.vd$labels <- c("CSF\n1333","nCSF\n1653","Other\n3834")
szpara.vd$colors <- c(0,0,0)

plot(szpara.vd)

