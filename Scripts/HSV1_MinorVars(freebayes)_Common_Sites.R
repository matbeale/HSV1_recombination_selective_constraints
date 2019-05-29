

# Use multiple sample outputs from MinorAlleleFreq_Freebayes.R

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

#basedropbox <- "C:/Users/Mat/Dropbox"
basedropbox <- "/Users/mb29/Dropbox"

#currentseq <- "HSV1-CSF1"
#seqlist <- c("HSV1-CSF1","HSV1-CSF2","HSV1-CSF3")
seqlist <- c("HSV1-CSF10","HSV1-CSF12","HSV1-CSF13","HSV1-CSF1","HSV1-CSF2","HSV1-CSF3","HSV1-CSF4","HSV1-CSF5","HSV1-nCSF10","HSV1-nCSF11","HSV1-nCSF13","HSV1-nCSF14","HSV1-nCSF15","HSV1-nCSF1","HSV1-nCSF4","HSV1-nCSF5","HSV1-nCSF6","HSV1-nCSF7")

minorvars.collating <- data.frame(c(1:152222),stringsAsFactors = F)
colnames(minorvars.collating) <- "position"

for (currentseq in seqlist) {
seqpath <- paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/")
#seqpath <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/"

seqname <- paste(currentseq,"-freebayes.MinorVars.freebayes.txt", sep="")
minorvars.called <- read.table(paste(seqpath,seqname,sep=""), header=T, stringsAsFactors = F)
minorvars.called <- data.frame(minorvars.called$position,stringsAsFactors = F)
minorvars.called$seq <- minorvars.called[,1]
colnames(minorvars.called)[ncol(minorvars.called)] <- currentseq


minorvars.collating <- merge(minorvars.collating,minorvars.called, by.x="position",by.y="minorvars.called.position",all.x=T)
#minorvars.collating[ncol(minorvars.collating)] <- sapply(1:nrow(minorvars.collating), function(x) ifelse(is.na(minorvars.collating[x,ncol(minorvars.collating)]),0,minorvars.collating[x,ncol(minorvars.collating)]))
minorvars.collating[ncol(minorvars.collating)] <- sapply(1:nrow(minorvars.collating), function(x) ifelse(is.na(minorvars.collating[x,ncol(minorvars.collating)]),0,1))
}

####
# Write to table and/or read in to avoid above loop (slow)
#write.table(minorvars.collating, file=paste(seqpath,"All-MinorVar-sites.freebayes.txt", sep=""), quote=T, sep="\t",col.names=T, row.names=F)
minorvars.collating <- read.table(paste(seqpath,"All-MinorVar-sites.freebayes.txt", sep=""), stringsAsFactors = F, header=T)
colnames(minorvars.collating) <- gsub("\\.","\\-",colnames(minorvars.collating))


# strip out positions without any minor variants in dataset
minorvars.collating$keep <- sapply(1:nrow(minorvars.collating), function(x) sum(minorvars.collating[x,2:ncol(minorvars.collating)]))
minorvars.collating <- subset(minorvars.collating, keep>0, select=-keep)


# strip out massively mixed sample (HSV1-nCSF7)
minorvars.collating5 <- minorvars.collating # keep all samples in case needed later
minorvars.collating<-minorvars.collating[,c(1:18)]


# Plot tile layout
minorvars.melted <- melt(minorvars.collating, id.vars="position")
minorvars.melted$source <- gsub("HSV1-","",minorvars.melted$variable)
minorvars.melted$source <- gsub("CSF.+$","CSF",minorvars.melted$source,perl=T)
minorvars.melted$minors <- ifelse(minorvars.melted$value >0, "minor","none")
minorvars.melted$position <- as.factor(minorvars.melted$position)

p2 <- ggplot(minorvars.melted) + facet_grid(source ~ ., scales="free", space="free", drop=T, margins=F)
p2 <- p2 + geom_tile(aes(x=position, y=variable, fill=minors)) + 
  #geom_tile(aes(x=genepos, y=syn, label=syn)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme(legend.position = "right") + labs(x="SNP Position in Genome", y="Sequence") +
  scale_fill_discrete(name="MinorVar") #+ ggtitle(currentgene)
p2 

Cairo(file=paste(seqpath,"All-MinorVar-sites.freebayes.tiling.png",sep=""), width = 1800, height = 300,type="png",dpi=1200, units = "pt")
#paste(seqpath,"All-MinorVar-sites.freebayes.tiling.png",sep="")
p2
dev.off()

# Look only at minor-var sites occurring in multiple genomes
minorvars.collating3 <- minorvars.collating
minorvars.collating3$multiple <- sapply(1:nrow(minorvars.collating3), function(x) sum(minorvars.collating3[x,2:ncol(minorvars.collating)]))
minorvars.collating3 <- subset(minorvars.collating3, multiple>1, select=-multiple)
# Plot tile layout
minorvars.melted2 <- melt(minorvars.collating3, id.vars="position")
minorvars.melted2$source <- gsub("HSV1-","",minorvars.melted2$variable)
minorvars.melted2$source <- gsub("CSF.+$","CSF",minorvars.melted2$source,perl=T)
minorvars.melted2$minors <- ifelse(minorvars.melted2$value >0, "minor","none")
minorvars.melted2$position <- as.factor(minorvars.melted2$position)
p3 <- ggplot(minorvars.melted2) + facet_grid(source ~ ., scales="free", space="free", drop=T, margins=F)
p3 <- p3 + geom_tile(aes(x=position, y=variable, fill=minors)) + 
  #geom_tile(aes(x=genepos, y=syn, label=syn)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme(legend.position = "right") + labs(x="SNP Position in Genome", y="Sequence") +
  scale_fill_discrete(name="MinorVar") #+ ggtitle(currentgene)
p3
Cairo(file=paste(seqpath,"All-MinorVar--multihit-sites.freebayes.tiling.png",sep=""), width = 1200, height = 300,type="png",dpi=1200, units = "pt")
p3
dev.off()
#??? plot this with upsetR ???

write.table(minorvars.collating3, file=paste(seqpath,"All-MinorVar-multihit-sites.freebayes.tiling.txt",sep=""),quote=T, sep="\t",col.names=T, row.names=F)





# pull in coding information & gene positions
gene.matrix <- read.table(paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1_gene-matrix.txt"),header =T, stringsAsFactors = F)
#gene.matrix <- read.table("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1_gene-matrix.txt",header =T, stringsAsFactors = F)
gene.matrix.melted <- melt(gene.matrix, id.vars="Position")
gene.matrix.melted$value <- as.factor(gene.matrix.melted$value)
# bring in gene names for text and determine midpoint in gene position for plotting
## Pull in gene positions for mapping hits
gene.pos.file <- paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified")
#gene.pos.file <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified"
gene.pos <- read.table(gene.pos.file, header =F, stringsAsFactors = F, quote="")
colnames(gene.pos) <- c("Start","End","Gene","Strand")
gene.pos$Start <- as.numeric(gene.pos$Start)
gene.pos$End <- as.numeric(gene.pos$End)
gene.pos$midpoint <- as.numeric(gene.pos$Start) + ((as.numeric(gene.pos$End) - as.numeric(gene.pos$Start) )/2)


######
# Assign each variant to gene
minorvars.collating2 <- minorvars.collating
minorvars.collating2$genemappings <- gsub("^, ","",sapply(1:nrow(minorvars.collating2), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(minorvars.collating2$position[y] >= gene.pos$Start[x] && minorvars.collating2$position[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)

# Try and look at coding (syn/non-syn)
### this is the old Core snps db### 
#coding.df <- read.table(file="C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/HSV1_Core-SNPs-dedup.coding-data", header=T, stringsAsFactors=F)
# this is the new one with all snps including minor var sites included###

coding.df <- read.table(file="C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/HSV1_MinorVars-freebayes-AllSNPs.coding-data", header=T, stringsAsFactors=F)#, col.names=T)#, row.names=F)

minorvars.collating2 <- merge(minorvars.collating2, coding.df, by.x="position", by.y="start", all.x=T)

# Now look at only the non-syn variants
minorvars.collating2$CONSEQUENCE <- sapply(1:nrow(minorvars.collating2), function(x) ifelse(is.na(minorvars.collating2$CONSEQUENCE[x]),"",minorvars.collating2$CONSEQUENCE[x]))
minorvars.nonsyn <- subset(minorvars.collating2, CONSEQUENCE=="nonsynonymous", select=c(1:19))

# Further filter to look at sites affecting multiple sequences
minorvars.nonsyn$multiple <- as.numeric(sapply(1:nrow(minorvars.nonsyn), function(x) sum(minorvars.nonsyn[x,2:(ncol(minorvars.nonsyn)-1)])))
minorvars.nonsyn <- subset(minorvars.nonsyn, multiple > 1, select=-multiple)

#write.table(minorvars.nonsyn, file=paste(seqpath,"All-MinorVar-nonSyn-multihit-sites.freebayes.tiling.txt",sep=""),quote=T, sep="\t",col.names=T, row.names=F)

minorvars.nonsyn <- read.table(file=paste(seqpath,"All-MinorVar-nonSyn-multihit-sites.freebayes.tiling.txt",sep=""), header=T, stringsAsFactors=F)#, col.names=T)#, row.names=F)
colnames(minorvars.nonsyn) <- gsub("\\.","\\-",colnames(minorvars.nonsyn))

# Plot tile layout
minorvars.melted3 <- melt(minorvars.nonsyn, id.vars="position")
minorvars.melted3$source <- gsub("HSV1-","",minorvars.melted3$variable)
minorvars.melted3$source <- gsub("CSF.+$","CSF",minorvars.melted3$source,perl=T)
minorvars.melted3$minors <- ifelse(minorvars.melted3$value >0, "minor","none")
minorvars.melted3$position <- as.factor(minorvars.melted3$position)
#p4 <- ggplot(minorvars.melted3) + facet_grid(source ~ ., scales="free", space="free", drop=T, margins=F)
#p4 <- p4 + geom_tile(aes(x=position, y=variable, fill=minors)) + 
  #geom_tile(aes(x=genepos, y=syn, label=syn)) + 
#  theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
#  theme(legend.position = "top") + labs(x="SNP Position in Genome", y="Sequence") +
#  scale_fill_discrete(name="MinorVar") #+ ggtitle(currentgene)
#p4

###
#try adding gene labels#
p4 <- ggplot(minorvars.melted3) 
p4 <- p4 + geom_tile(data=subset(minorvars.melted3,variable != "genemappings"), aes(x=position, y=variable, fill=minors)) + 
  #geom_tile(aes(x=genepos, y=syn, label=syn)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme(legend.position = "top") + labs(x="SNP Position in Genome", y="Sequence") +
  scale_fill_discrete(name="MinorVar") #+ scale_color_manual(values=c("red","white") )
p4 <- p4 + geom_text(data=subset(minorvars.melted3,variable == "genemappings"), aes(x=position, y="Gene", label=value), angle = -90, size=2.5, alpha=1/3) + facet_grid(source ~ ., scales="free", space="free", drop=T, margins=F)

Cairo(file=paste(seqpath,"All-MinorVar-nonSyn-sites.freebayes.tiling.png",sep=""), width = 1200, height = 700,type="png",dpi=1200, units = "pt")
#Cairo(file=paste(seqpath,"All-MinorVar-multihit-nonSyn-sites.freebayes.tiling.png",sep=""), width = 1200, height = 400,type="png",dpi=1200, units = "pt")
p4
dev.off()




#########################
################
# Evaluate all minor variant sites based on coverage at that position in all samples
highcov.file <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_250x_0.9_minorvar.freebayes.fullcount"
highcov.names <- read.table(highcov.file, header=F, sep="")
highcov.names <- as.character(unlist(highcov.names[,1]))
write.table(highcov.file, paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_250x_0.9_minorvar.freebayes.fullcount"))

# Samples with >250x at >90% genome
#highcov.names <- c("HSV1-nCSF11","HSV1-nCSF10","HSV1-nCSF7","HSV1-nCSF14","HSV1-nCSF15","HSV1-nCSF5","HSV1-CSF1","HSV1-nCSF1","HSV1-CSF4","HSV1-CSF2","HSV1-CSF5")
# Samples with >250x at >85% genome
highcov.names <- c("HSV1-nCSF11","HSV1-nCSF10","HSV1-nCSF7","HSV1-nCSF14","HSV1-nCSF15","HSV1-nCSF5","HSV1-CSF1","HSV1-nCSF1","HSV1-CSF4","HSV1-CSF2","HSV1-CSF5","HSV1-nCSF6","HSV1-CSF3","HSV1-nCSF13")

# Samples with >250x at >85% genome including duplicates
#highcov.names <- c("HSV1-nonCSF13",	"HSV1-CSF3",	"HSV1-nonCSF6",	"HSV1-CSF2",	"HSV1-CSF5",	"HSV1-CSF4",	"HSV1-nonCSF1",	"HSV1-nonCSF3",	"HSV1-CSF1",	"HSV1-nonCSF5",	"HSV1-nonCSF9",	"HSV1-nonCSF14",	"HSV1-nonCSF15",	"HSV1-nonCSF7",	"HSV1-nonCSF10",	"HSV1-nonCSF11")



sample.threshold <- length(highcov.names)

minorvars.cov <- data.frame(minorvars.collating$position, stringsAsFactors = F)
colnames(minorvars.cov) <- "position"


####
# Need to make a copy of either the raw coverage into Dropbox OR the output of this loop

for (currentseq in highcov.names) {  
  covpath <- "C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/Coverage/RawCoverage/"
  seqname <- paste(currentseq,".depth", sep="")
  coverage.called <- read.table(paste(covpath,seqname,sep=""), header=T, stringsAsFactors = F)
  colnames(coverage.called) <- c("ref","position",paste(currentseq,".cov",sep=""))
  coverage.called <- coverage.called[,-1]
  # collate coverage at each minor site
  minorvars.cov <- merge(minorvars.cov, coverage.called, by="position")
  # evaluate if site meets minimum threshold (250x)
  minorvars.cov[ncol(minorvars.cov)] <- sapply(1:nrow(minorvars.cov), function(x) ifelse(minorvars.cov[x,ncol(minorvars.cov)] <250,0,1))
}
minorvars.cov$passed <- sapply(1:nrow(minorvars.cov), function(x) sum(minorvars.cov[x,c(2:ncol(minorvars.cov))]))


#write.table(minorvars.cov, paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_Minorvars.Highcov.txt"))
minorvars.cov <- read.table(paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_Minorvars.Highcov.txt"))

#####

# Only use highcov
highcov.minorvars <- minorvars.collating5[,colnames(minorvars.collating5)%in%highcov.names]

highcov.minorvars$position <- minorvars.collating5$position
minorvars.cov.passed <- data.frame(minorvars.cov[,c("position","passed")],stringsAsFactors = F)
highcov.minorvars <- merge(highcov.minorvars, minorvars.cov.passed, by="position")
nrow(subset(highcov.minorvars, passed==sample.threshold, select=-passed))

highcov.minorvars <- subset(highcov.minorvars, passed==sample.threshold, select=-passed)

#highcov.minorvars.counts <- data.frame(highcov.names)
highcov.minorvars.counts <- NULL
for (sequence in highcov.names){
  sum(highcov.minorvars[,sequence])
  highcov.minorvars.counts <- cbind(highcov.minorvars.counts,sum(highcov.minorvars[,sequence]))
  colnames(highcov.minorvars.counts)[ncol(highcov.minorvars.counts)] <- sequence
}
highcov.minorvars.counts <- data.frame(t(highcov.minorvars.counts),stringsAsFactors = F)
highcov.minorvars.counts$Sample <- rownames(highcov.minorvars.counts)
colnames(highcov.minorvars.counts) <- c("freq","Sample")
highcov.minorvars.counts$SampleSite <- ifelse(grepl("nCSF",highcov.minorvars.counts$Sample),"nCSF","CSF")                                    
highcov.minorvars.counts$order <- xtfrm(highcov.minorvars.counts$freq)
highcov.minorvars.counts <- highcov.minorvars.counts[with(highcov.minorvars.counts, order(order, Sample)), ]
highcov.minorvars.counts$order3 <- factor(highcov.minorvars.counts$Sample, as.character(highcov.minorvars.counts$Sample))

# Now plot
highcov.minorvars.counts2 <- highcov.minorvars.counts
highcov.minorvars.counts2$SampleSite <- gsub("nCSF","SWAB",highcov.minorvars.counts2$SampleSite)

#p19 <- ggplot(highcov.minorvars.counts, aes(x=order3,y=freq, fill=SampleSite))
p19 <- ggplot(highcov.minorvars.counts2, aes(x=order3,y=freq, fill=SampleSite))
p19 <- p19 + geom_bar(stat="identity") +  
  theme_bw() + 
  labs(x="Sequence", y="Minor Variant Sites\n(count/genome)") +
  #scale_y_log10(expand=c(0,0), breaks=c(1,5,10,20, 50,100, max(highcov.minorvars.counts$freq)), limits=c(1,200)) +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(0, 150,25)), limits=c(0,160)) +
  #scale_fill_manual(values = c("#009900","red3")) +
  scale_fill_manual(values = c("#009900","darkorchid3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5), legend.position="top")
p19

# Use a boxplot
highcov.minorvars.counts$freq <- as.numeric(highcov.minorvars.counts$freq)
p17 <- ggplot(highcov.minorvars.counts, aes(SampleSite, freq,colour=SampleSite))
#p17 <- p17 + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","red3")) +
p17 <- p17 + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","darkorchid3")) +  
  labs(x="Sample Site", y="Minor Variant Sites (count/genome)") +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(0, 150,25)), limits=c(0,160))
  #scale_y_log10(expand=c(0,0), breaks=c(1,5,10,20, 50,100, max(highcov.minorvars.counts$freq)), limits=c(1,200))
p17


# Combine sample plot with box plot
gA <- ggplotGrob(p19 + theme(legend.justification=c(0,0), legend.position=c(0.1,0.75) ,legend.text = element_text(size = 10 ),legend.key.size = unit(0.5,"cm"),legend.direction="horizontal",axis.title.x = element_blank() ))
gB <- ggplotGrob(p17 + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank(),axis.title.x = element_blank() ))
gB$heights <- gA$heights
grid.arrange(gA,gB, nrow = 1, widths=c(4,1))

#Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_250x_0.85(364sites)_minorvar.freebayes.fullcount.png", width = 800, height = 300,type="png",dpi=900, units = "pt")
Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_250x_0.85(364sites)_minorvar_nolog.freebayes.fullcount.png", width = 800, height = 300,type="png",dpi=900, units = "pt")
grid.arrange(gA,gB, nrow = 1, widths=c(4,1))
dev.off()

##############################################

# Repeat Minor Var analysis, but keep frequency data
minorvars.collating.freq <- data.frame(c(1:152222),stringsAsFactors = F)
colnames(minorvars.collating.freq) <- "position"

for (currentseq in seqlist) {
  seqpath <- paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/")
  seqname <- paste(currentseq,"-freebayes.MinorVars.freebayes.txt", sep="")
  minorvars.called <- read.table(paste(seqpath,seqname,sep=""), header=T, stringsAsFactors = F)
  minorvars.called <- data.frame(minorvars.called[c("position","minorfreq")],stringsAsFactors = F)
  colnames(minorvars.called)[ncol(minorvars.called)-1] <- "position"
  colnames(minorvars.called)[ncol(minorvars.called)] <- currentseq #paste(currentseq,".freq", sep="")
  
  minorvars.collating.freq <- merge(minorvars.collating.freq,minorvars.called, by="position",all.x=T)
  minorvars.collating.freq[ncol(minorvars.collating.freq)] <- sapply(1:nrow(minorvars.collating.freq), function(x) ifelse(is.na(minorvars.collating.freq[x,ncol(minorvars.collating.freq)]),0,minorvars.collating.freq[x,ncol(minorvars.collating.freq)]))
  #minorvars.collating[ncol(minorvars.collating)] <- sapply(1:nrow(minorvars.collating), function(x) ifelse(is.na(minorvars.collating[x,ncol(minorvars.collating)]),0,1))
}

####
# Write to table and/or read in to avoid above loop (slow)
#write.table(minorvars.collating.freq, file=paste(seqpath,"All-MinorVar-freqs.freebayes.txt", sep=""), quote=T, sep="\t",col.names=T, row.names=F)
minorvars.collating.freq <- read.table(paste(seqpath,"All-MinorVar-freqs.freebayes.txt", sep=""), stringsAsFactors = F, header=T)
colnames(minorvars.collating.freq) <- gsub("\\.","\\-",colnames(minorvars.collating.freq))

# strip out positions without any minor variants in dataset
minorvars.collating.freq$keep <- sapply(1:nrow(minorvars.collating.freq), function(x) sum(minorvars.collating.freq[x,2:ncol(minorvars.collating.freq)]))
minorvars.collating.freq <- subset(minorvars.collating.freq, keep>0, select=-keep)

# strip out massively mixed sample (HSV1-nCSF7)
minorvars.collating.freq2 <- minorvars.collating.freq # keep all samples in case needed later
#minorvars.collating.freq<-minorvars.collating.freq[,c(1:18)]

# only include sites previously determined to be covered at >250x in all samples
minorvars.collating.freq2 <- minorvars.collating.freq2[,colnames(minorvars.collating.freq2)%in%highcov.names]
minorvars.collating.freq2 <- minorvars.collating.freq2[rownames(minorvars.collating.freq2)%in%highcov.minorvars$position,]

# Generate freq bins and bin variants
binsize <- 2
minorvars.freq.bins <- data.frame(rownames(minorvars.collating.freq2),stringsAsFactors = F)
colnames(minorvars.freq.bins) <- "position"
for (currentseq in colnames(minorvars.collating.freq2)){
  #minorvars.freq.bins$freqbin <- ((ceiling(minorvars.collating.freq2[which(colnames(minorvars.collating.freq2)==currentseq)] / binsize))*binsize)-(binsize/2) # create binning variables and re-label to nearest bin 
  minorvars.freq.bins$freqbin <- as.numeric(unlist(((ceiling(minorvars.collating.freq2[which(colnames(minorvars.collating.freq2)==currentseq)] / binsize))*binsize)-(binsize/2)))
  colnames(minorvars.freq.bins)[ncol(minorvars.freq.bins)] <- currentseq
}

# now do a melt and assign sample source for colouring
minorvars.freq.bins.melted <- melt(minorvars.freq.bins,id.var="position")
minorvars.freq.bins.melted$Source <- ifelse(grepl("nCSF",minorvars.freq.bins.melted$variable),"nCSF","CSF")

### Reorder data based on those used in total freq plot
minorvars.freq.bins.melted$order <- minorvars.freq.bins.melted$variable
minorvars.freq.bins.melted$order <- gsub("nCSF","SWAB", minorvars.freq.bins.melted$order)

#minorvars.freq.bins.melted$order<- factor(minorvars.freq.bins.melted$order,levels=as.character(unlist(highcov.minorvars.counts$order3)))

minorvars.freq.bins.melted$order<- factor(minorvars.freq.bins.melted$order,levels=gsub("nCSF","SWAB",as.character(unlist(highcov.minorvars.counts$order3))))
#gsub("nCSF","SWAB",as.character(unlist(highcov.minorvars.counts$order3)))

#####

# Change sample names (nCSF to SWAB)



# now plot with facets...
p3 <- ggplot(minorvars.freq.bins.melted, aes(value)) #+ facet_grid(variable)
p3 <- p3 + geom_bar(stat="count",aes(fill=Source)) +
#p3 <- p3 + geom_jitter(stat="count",aes(colour=Source)) +  
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  labs(x="", y="Occurrences in Genome (count)") +
  scale_y_continuous(expand=c(0,0), limits=c(0,60)) +
  #scale_y_log10(expand=c(0,0), limits=c(1,50)) +
  scale_x_continuous( expand=c(0,0),breaks=seq(0,50,binsize),limits = c(0, 50)) + # #need to specify binsize in limits or 48-50 is excluded
  theme(axis.text.x = element_text(hjust = 1, vjust=0.5,angle=90),strip.text.x = element_text(angle=90)) +
  coord_flip() + scale_fill_manual(values = c("#009900","red3")) +
  #facet_wrap(~ variable)#, scales="free", space="free", drop=T, margins=F)
  facet_grid(. ~ order, scales="free", space="free", drop=T, margins=F)
  #theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())
p3


#Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar.freebayes.histogram-plots.png", width = 1400, height = 300,type="png",dpi=1200, units = "pt")
#Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.histogram-plots.png", width = 1400, height = 300,type="png",dpi=1200, units = "pt")
Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites)-bin2.freebayes.histogram-plots.png", width = 1400, height = 300,type="png",dpi=1200, units = "pt")

p3
dev.off()

# Density plot
p3 <- ggplot(minorvars.freq.bins.melted, aes(value)) #+ facet_grid(variable)
p3 <- p3 + geom_density(stat="density",aes(alpha=1/50, group=variable,colour=variable, fill=variable),position = "stack") +
  #p3 <- p3 + geom_jitter(stat="count",aes(colour=Source)) +  
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw() +
  labs(x="Minor Variant Frequency", y="Occurrences in Genome (count)") +
  scale_x_continuous(limits = c(1, 50)) #+ # #need to specify binsize in limits or 48-50 is excluded

p3
Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.density-plots.png", width = 1400, height = 300,type="png",dpi=1200, units = "pt")
p3
dev.off()


# Violin Plot
p3 <- ggplot(minorvars.freq.bins.melted, aes(factor(order),value)) #+ facet_grid(variable)
p3 <- p3 + geom_violin(scale = "area", adjust = 0.5,aes(alpha=1/20, group=variable, fill=Source, colour=Source)) +
  theme_bw() + ylim(0,50) + 
  scale_fill_manual(values = c("#009900","red3")) + scale_colour_manual(values = c("#009900","red3")) +
  geom_jitter(height=0, alpha=1/2) +  
  #geom_bar(stat="count",alpha=1/2) +  
  theme(axis.text.x = element_text(hjust = 1, vjust=0.5,angle=-90)) +
  labs(y="Minor Variant Frequency & Density", x="Sample")
p3  
Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.violin-plots.png", width = 900, height = 300,type="png",dpi=900, units = "pt")
p3
dev.off()  




# Use geom_sina for plot ???
library("ggforce")



# Change sample names (nCSF to SWAB)
#minorvars.freq.bins.melted$order <- gsub("nCSF","SWAB", minorvars.freq.bins.melted$order)



p3 <- ggplot(minorvars.freq.bins.melted, aes(factor(order),value)) #+ facet_grid(variable)
#p3 <- p3 + geom_sina( alpha=1/2, size=2, shape=3,aes( group=variable, fill=Source, colour=Source)) +

#p3 <- p3 + geom_sina( alpha=1/2, size=2, aes( group=variable, fill=Source, colour=Source)) +  
p3 <- p3 + geom_sina( alpha=1/4, size=2, binwidth=5,scale=F,method ="counts", maxwidth=1, aes( group=variable, fill=Source, colour=Source)) +  
  theme_bw() + ylim(0,50) + #coord_cartesian(ylim = c(0, 50)) + 
  #scale_fill_manual(values = c("#009900","red3")) + scale_colour_manual(values = c("#009900","red3")) +
  scale_fill_manual(values = c("#009900","darkorchid3")) + scale_colour_manual(values = c("#009900","darkorchid3")) +
  #geom_jitter(height=0, alpha=1/2) +  
  #geom_bar(stat="count",alpha=1/2) +  
  theme(axis.text.x = element_text(hjust = 1, vjust=0.5,angle=90)) +
  #labs(y="Minor Variant\nFrequency & Density", x="Sample")# +
  labs(y="Variant Frequency (%)", x="Sample")# +
  #facet_zoom(x = order == "HSV1-nCSF7")
p3  
Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.Sina-plots.png", width = 900, height = 300,type="png",dpi=900, units = "pt")
p3
dev.off() 







# Combine sample plot with box plot
gA <- ggplotGrob(p19 + theme(legend.justification=c(0,0), legend.position=c(0.1,0.75) ,legend.text = element_text(size = 10 ),legend.key.size = unit(0.5,"cm"),legend.direction="horizontal",axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank()))
gB <- ggplotGrob(p17 + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank(),axis.title.x = element_blank() ))
gC <- ggplotGrob(p3 + theme(legend.position="none"))
gB$heights <- gA$heights
gC$widths <- gA$widths
grid.arrange(gA,gB,gC,nrow = 2, widths=c(4,1), heights=c(1,2))

Cairo(file= "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.Sina-plots+Counts.png", width = 450, height = 300,type="png",dpi=600, units = "pt")
grid.arrange(gA,gB,gC,nrow = 2, widths=c(4,1), heights=c(2,3))
dev.off() 




# attempt to grid plots together - not yet working, as samples need reordering
#grid.arrange(p19,p3, nrow=2)



#####
# Pull in Juli Cudini's analysis of intrahost diversity (nuc div)
juli.file <- paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/Juli_Cudini_Intrahost_diversity/","HSV_bootstap.matpaper.tsv")
#juli.nucdiv <- read.csv(juli.file, stringsAsFactors = F, header = T)
juli.nucdiv <- read.table(juli.file, stringsAsFactors = F, header = T)

colnames(juli.nucdiv) <- gsub("\\.","-",colnames(juli.nucdiv))
colnames(juli.nucdiv) <- gsub("non","n",colnames(juli.nucdiv))
colnames(juli.nucdiv) <- gsub("HSV-","HSV1-",colnames(juli.nucdiv))

#juli.nucdiv.sub <- juli.nucdiv[,c(1,highcov.names)]
juli.nucdiv.sub <- juli.nucdiv[,c("SIM",highcov.names)]
juli.melt <- melt(juli.nucdiv, id.vars="SIM")

# reorder
juli.melt$order <- juli.melt$variable
juli.melt$order<- factor(juli.melt$order,levels=as.character(unlist(highcov.minorvars.counts$order3)))
#
# Group by CSF/nonCSF
juli.melt$Source <- gsub("HSV1-","",gsub("CSF.+$","CSF",juli.melt$variable))
juli.melt$Source <- gsub("nCSF","SWAB",juli.melt$Source)


p.juli <- ggplot(juli.melt, aes(factor(order), value, colour=Source)) 
p.juli <- p.juli + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","darkorchid3")) +
  #scale_colour_manual(values = c("springgreen4","red4")) +
  labs(x="Sample", y=expression(paste("Nucleotide Diversity (", pi, ")")))
p.juli

####
# Box plot by sample site
p.juli.group <- ggplot(juli.melt, aes(Source, value, colour=Source)) 
#p.juli.group <- p.juli.group + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","red3")) +
p.juli.group <- p.juli.group + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","darkorchid3"))  
  #scale_colour_manual(values = c("springgreen4","red4")) +
  labs(x="Sample Site", y=expression(paste("Intra-host nucleotide diversity (", pi, ")")))
p.juli.group


# Combine sample plot with box plot
gA <- ggplotGrob(p19 + theme(legend.justification=c(0,0), legend.position=c(0.1,0.75) ,legend.text = element_text(size = 10 ),legend.key.size = unit(0.5,"cm"),legend.direction="horizontal",axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank()))
gB <- ggplotGrob(p17 + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank(), axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank() ))
gC <- ggplotGrob(p.juli + theme(legend.position="none",axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank()))
#gD <- ggplotGrob(p.juli.group + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank()))
gD <- ggplotGrob(p.juli.group + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank()))
gE <- ggplotGrob(p3 + theme(legend.position="none"))

#gr <- ggplotGrob(r)
gB$heights <- gA$heights
#gC$heights <- gD$heights
gD$heights <- gC$heights
gC$widths <- gA$widths
gD$widths <- gB$widths
gE$widths <- gA$widths

#r$heights <- gC$heights
grid.arrange(gA,gB,gC,gD,gE,ncol=2,nrow = 3, widths=c(4,1), heights=c(1,1,2))


Cairo(file= paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.Sina-plots+Counts+JuliPi.png"), width = 500, height = 450,type="png",dpi=600, units = "pt")
grid.arrange(gA,gB,gC,gD,gE,ncol=2,nrow = 3, widths=c(4,1), heights=c(2,2,3))
dev.off() 







# what about a boxplot for minorvar frequency aswell?
minorvars.freq.bins.melted.varonly <- minorvars.freq.bins.melted[minorvars.freq.bins.melted$value>=2,]
p3.source <- ggplot(minorvars.freq.bins.melted.varonly, aes(Source,value,colour=Source)) 
p3.source <- p3.source + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","red3")) +
  #scale_colour_manual(values = c("springgreen4","red4")) +
  labs(x="Sample Site", y="Minor Variant\nFrequency & Density")
p3.source





# Combine sample plot with box plot
gA <- ggplotGrob(p19 + theme(legend.justification=c(0,0), legend.position=c(0.1,0.75) ,legend.text = element_text(size = 10 ),legend.key.size = unit(0.5,"cm"),legend.direction="horizontal",axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank()))
gB <- ggplotGrob(p17 + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank(), axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank() ))
gC <- ggplotGrob(p.juli + theme(legend.position="none",axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank()))
gD <- ggplotGrob(p.juli.group + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank(),axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x= element_blank() ))
gE <- ggplotGrob(p3 + theme(legend.position="none"))
gF <- ggplotGrob(p3.source + theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank()))

#gr <- ggplotGrob(r)
gB$heights <- gA$heights
gC$heights <- gD$heights
gF$heights <- gE$heights
gC$widths <- gA$widths
gD$widths <- gB$widths
gE$widths <- gA$widths
gF$widths <- gB$widths

#r$heights <- gC$heights
grid.arrange(gA,gB,gC,gD,gE,gF,ncol=2,nrow = 3, widths=c(4,1), heights=c(1,1,2))


Cairo(file= paste0(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/MinorVars_(ref-based)/FreeBayes_calling/HSV1_all-minorvar_250x(364sites).freebayes.Sina-plots+Counts+JuliPi+freqdensity.png"), width = 450, height = 450,type="png",dpi=600, units = "pt")
grid.arrange(gA,gB,gC,gD,gE,gF,ncol=2,nrow = 3, widths=c(4,1), heights=c(1,1,2))
dev.off() 






# Look at p3.source again - is it biased by nCSF7 (loads of vars throughout) - remove it and replot
# Could also be biased by CSF3 (loads of low freq vars)
# what about a boxplot for minorvar frequency aswell?

minorvars.freq.bins.melted.varonly.rnCSF7 <- minorvars.freq.bins.melted.varonly[minorvars.freq.bins.melted.varonly$variable!="HSV1-nCSF7",]
minorvars.freq.bins.melted.varonly.rnCSF7 <- minorvars.freq.bins.melted.varonly.rnCSF7[minorvars.freq.bins.melted.varonly.rnCSF7$variable!="HSV1-CSF3",]

p3.source.rnCSF7 <- ggplot(minorvars.freq.bins.melted.varonly.rnCSF7, aes(Source,value,colour=Source)) 
p3.source.rnCSF7 <- p3.source.rnCSF7 + geom_boxplot() + theme_bw() +scale_colour_manual(values = c("#009900","red3")) +
  #scale_colour_manual(values = c("springgreen4","red4")) +
  labs(x="Sample Site", y="Minor Variant\nFrequency & Density")
p3.source.rnCSF7
# If we take out nCSF7 and CSF3, the pattern goes away. I think CSF3 is the major bias here. 
