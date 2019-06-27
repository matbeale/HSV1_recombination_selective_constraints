#!/usr/bin/env Rscript

# This script takes in a multiple sequence alignment in .fasta format and generates plots
# that show the position of any variants (compared to the first genome in the file, 
# e.g. a reference). Since genomes of over ~5Kb will have problems with a 'plot all sites' 
# approach because there are too many sites to fit into a single image, I have included 
# a second plot that implements sliding windows to show whether a SNP occurs within a 
# defined window (e.g. 10bp).  


# Code to run from the command line is deprecated and commented out. To run in Rstudion the 
# input fasta file can be specified using the 'genefile' variable.
################
# Usage : GW_snps_plot.R \
#   <multiple-seq-alignment.fas>           # 1

#args <- commandArgs(trailingOnly = TRUE) # capture command line arguments
# Put command line arguments into variables
#genefile <- args[1]
###################

# Manually specify path to fasta file input

#genefile <- "/Users/mb29/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Gubbins/Szpara2014-v4/HSV1_Corededup+Szpara2014_v3_.fas"
genefile <- "/Users/mb29/hsv1/parsnp/assembly-vs-Szpara/parsnp/all_genomes/iqtree/HSV1-assemblies-vs-Szpara.mafft-add.fa"

#genefile <- "/Users/mb29/Dropbox/UCL_Pathseek/HBV/HBV-Genes/HBVall_S.Snt.fasta"
#genefile <- "/Users/mb29/Dropbox/UCL_Pathseek/HSV1/Ref_genomes/HSV1_vs_HSV2/HSV1-vs-HSV2_ref.aln.fas"
#genefile <- "/Users/mb29/Dropbox/UCL_Pathseek/HSV1/Ref_based/raxml/HSV1-snippy-core/CSF_possible-same-patient.fas"

#genefile <- "/Users/mb29/Treponema/Lukehart_UoW/MiniProjects1_intrahost/Analysis/Phylo/SLK3a_bwa.aln"

##################


# Load Dependencies
if(!require(seqinr)){
  install.packages("seqinr",repos="http://www.stats.bris.ac.uk/R/")
  library("seqinr")
}
if(!require(ggplot2)){
  install.packages("ggplot2",repos="http://www.stats.bris.ac.uk/R/")
  library("ggplot2")
}
if(!require(Cairo)){
  install.packages("Cairo",repos="http://www.stats.bris.ac.uk/R/")
  library("Cairo")
}
if(!require(reshape2)){
  install.packages("reshape2",repos="http://www.stats.bris.ac.uk/R/")
  library("reshape2")
}
if(!require(gridExtra)){
  install.packages("gridExtra",repos="http://www.stats.bris.ac.uk/R/")
  library("gridExtra")
}
if(!require(grid)){
  install.packages("grid",repos="http://www.stats.bris.ac.uk/R/")
  library("grid")
}
if(!require(ggExtra)){
  install.packages("ggExtra",repos="http://www.stats.bris.ac.uk/R/")
  library("ggExtra")
}

#######################
# Read in sequence alignment 
geneseqs <- read.fasta(file=genefile, seqtype="DNA", as.string=T,set.attributes=F)

# For capturing and relabelling the sequence names to make them pretty (this needs to be customised)
seqnames <- names(geneseqs)
seqnames <- gsub(".consensus.bak","",gsub(".mspades.herpesvir.binned","",gsub(".iva.herpesvir.binned","",gsub(".fa","",seqnames))))

# Convert sequence alignment file into a dataframe (gene.aln)
gene.aln <- data.frame(t(as.matrix.alignment(as.alignment(nb=length(seqnames),nam=NULL,seq=geneseqs,com=NULL))), stringsAsFactors = F)
colnames(gene.aln) <- seqnames

# Only keep samples of interest for HSV1 study
#gene.aln.study <- gene.aln[,c("Reference_Strain-17",colnames(gene.aln)[grepl("CSF",colnames(gene.aln))])]
gene.aln.study <- gene.aln[,c("NC001806.2_HSV1_s_17",colnames(gene.aln)[grepl("CSF",colnames(gene.aln))])]

# Relabel "nCSF" as "SWAB"
colnames(gene.aln.study) <- gsub("nCSF","SWAB",colnames(gene.aln.study))
gene.aln <- gene.aln.study

# Remove regions that are only gaps
gene.aln.nogap <- gene.aln[,c(colnames(gene.aln)[!grepl("NC001806.2_HSV1_s_17",colnames(gene.aln))])]
gapsites <- data.frame(position=c(1:nrow(gene.aln.nogap)),stringsAsFactors = F)
gapsites$alleles <- sapply(1:nrow(gapsites), function (x) unique(as.character(gene.aln.nogap[x,])))
gapsites.gaps <- gapsites[gapsites$alleles=="-",]
gene.aln.nogap <- gene.aln[gapsites[gapsites$alleles!="-","position"],]
  

# Identify windows of gaps

#library(dplyr)
gapsites.windows <- gapsites %>%
  dplyr::mutate(alleles = ifelse(alleles=="-", 1, 0)) %>%
  dplyr::group_by(group = cumsum(c(0, diff(alleles) != 0))) %>%
  dplyr::filter(alleles == 1 & n() > 1) %>%
  dplyr::summarize("Start.Site"=min(position),
            "End.Site"=max(position),
            "Length.of.Run"=n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(-matches("group"))
gapsites.windows <- data.frame(gapsites.windows,stringsAsFactors = F)


# Make a binary matrix from the gene alignment
gene.aln.binary <- numeric(nrow(gene.aln.nogap))
gene.aln.binary <- data.frame(t(rbind(gene.aln.binary,sapply(1:nrow(gene.aln.nogap), function(x) sapply(2:ncol(gene.aln.nogap), function(y) ifelse(gene.aln.nogap[x,y]==gene.aln.nogap[x,1],"0","1"))))),stringsAsFactors = F)
colnames(gene.aln.binary) <- colnames(gene.aln)

gene.aln.binary$position <- row.names(gene.aln.nogap)

gene.aln.binary.melt <- melt(gene.aln.binary, id.vars="position")
colnames(gene.aln.binary.melt) <- c("Position","Sequence","Variant")
gene.aln.binary.melt$Position <- as.numeric(gene.aln.binary.melt$Position)

# Generate plot
p1 <- ggplot(gene.aln.binary.melt) 
p1 <- p1 + geom_tile(aes(x=Position, y=Sequence, fill=Variant, color=NULL)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme(legend.position = "right") + labs(x="Genome Position", y="Sequence") +
  scale_color_manual(breaks=c("Reference", "Variant"),values=c("white", "#F8766D")) +
  scale_fill_manual(breaks=c("Reference", "Variant"),values=c("white", "#F8766D"))
p1  
#p1 +coord_cartesian(xlim=c(1800:1900)) # for subsetting a particular part of the plot


####################
# Aggregate SNPs into sliding windows for larger genomes
###
windowsize <- 1000       # specify size of window
genomelength <- nrow(gene.aln)
nwindows <- genomelength/windowsize

# Create a copy of binary matrix and define bins
gene.aln.binary2 <- gene.aln.binary
row.names(gene.aln.binary2) <- gene.aln.binary2$position
gene.aln.binary2$window <- ((trunc(as.numeric(rownames(gene.aln.binary2)) / windowsize,0))*windowsize) # create binning variables and re-label to nearest bin 

# Loop through all sequences, and if there is a SNP (1) within the window, make it a 1
seqnames2 <- colnames(gene.aln.binary2)[c(1:19)]
seqnames <- seqnames2
gene.aln.windowed <- unique(gene.aln.binary2$window)
snps.windowed.counts <- unique(gene.aln.binary2$window)

for (seq in seqnames2){
snps.window <- aggregate(as.numeric(gene.aln.binary2[,seq]), list(gene.aln.binary2$window), max) # aggregate depth based on bin and calculate mean
snps.window[,2] <- as.character(snps.window[,2])
gene.aln.windowed <- cbind(gene.aln.windowed,snps.window[,2])

snps.window.count <- aggregate(as.numeric(gene.aln.binary2[,seq]), list(gene.aln.binary2$window), sum) # aggregate depth based on bin and calculate mean
snps.window.count[,2] <- as.character(snps.window.count[,2])
snps.windowed.counts <- cbind(snps.windowed.counts,snps.window.count[,2])
}
gene.aln.windowed <- data.frame(gene.aln.windowed,stringsAsFactors = F)
colnames(gene.aln.windowed) <- c("window",seqnames)

snps.windowed.counts <- data.frame(snps.windowed.counts,stringsAsFactors = F)
colnames(snps.windowed.counts) <- c("window",seqnames)


# Generate the plot of SNP presence/absence
gene.aln.windowed.melt <- melt(gene.aln.windowed, id.vars="window")
colnames(gene.aln.windowed.melt) <- c("Position","Sequence","Variant")
gene.aln.windowed.melt$Position <- as.numeric(gene.aln.windowed.melt$Position)

p2 <- ggplot(gene.aln.windowed.melt) 
p2 <- p2 + geom_tile(aes(x=Position, y=Sequence, fill=Variant, color=NULL)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme(legend.position = "right") + labs(x="Genome Position", y="Sequence") +
  #scale_color_manual(breaks=c("Reference", "Variant"),values=c("white", "#F8766D")) +
  scale_fill_manual(breaks=c("Reference", "Variant"),values=c("white", "#F8766D")) +
  scale_x_continuous(labels = scales::comma,expand = c(0, 0), breaks=seq(0, genomelength, round(genomelength / 15,1 -nchar(ceiling(genomelength / 15))))) 
p2  
#p2 +coord_cartesian(xlim=c(1800:1900)) # for subsetting a particular part of the plot



####
# Generate a sliding window plot of SNP density (faceted by genome)
snps.windowed.counts.melt <- melt(snps.windowed.counts, id.vars="window")
colnames(snps.windowed.counts.melt) <- c("Position","Sequence","Count")
snps.windowed.counts.melt$Position <- as.numeric(snps.windowed.counts.melt$Position)
# Remove reference (since there are no SNPs of course)
snps.windowed.counts.melt <- snps.windowed.counts.melt[snps.windowed.counts.melt$Sequence!="NC001806.2_HSV1_s_17",]
# Bring in info about whether it's SWAB or CSF
snps.windowed.counts.melt$type <- ifelse(grepl("SWAB",snps.windowed.counts.melt$Sequence),"SWAB","CSF")


p.SNPdensity <- ggplot(snps.windowed.counts.melt,aes(Position,as.numeric(Count),fill=type)) + 
  geom_bar(stat="identity") +
  #facet_grid(Sequence~., scales="free") +
  facet_grid(Sequence~.) + 
  #coord_cartesian(ylim=c(0,100)) +
  coord_cartesian(xlim=c(0,genomelength)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  theme(legend.position = "right") + 
  labs(x=paste0("Genome Position"," (",windowsize," bp windows)"), y="SNP count",fill="Sample Site") +
  theme(strip.background = element_blank(),strip.text.y = element_text(angle=0)) + 
  scale_fill_manual(breaks=c("CSF", "SWAB"),values=c("#009900","darkorchid3")) +
  scale_x_continuous(labels = scales::comma,expand = c(0, 0), breaks=seq(0, genomelength, round(genomelength / 15,1 -nchar(ceiling(genomelength / 15))))) +
  #scale_y_continuous(expand = c(0, 0), breaks=c(0,50,100)) +
  NULL
p.SNPdensity
#p2 +coord_cartesian(xlim=c(1800:1900)) # for subsetting a particular part of the plot



# Show regions of no/poor coverage (gapsites)
p.SNPdensity <- p.SNPdensity + geom_rect(data=gapsites.windows, aes(xmin=Start.Site, xmax=End.Site, ymin=0, ymax=25), fill="blue", alpha=0.1,inherit.aes=F)
p.SNPdensity


# Bring in gene positions to plot
gene.pos.file <- paste("/Users/mb29/Dropbox","/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified",sep="")
gene.pos <- read.table(gene.pos.file, header =F, stringsAsFactors = F, quote="")
colnames(gene.pos) <- c("Start","End","Gene","Strand")
gene.pos$Start <- as.numeric(gene.pos$Start)
gene.pos$End <- as.numeric(gene.pos$End)
gene.pos$Start.dir <- ifelse(gene.pos$Strand=='+', gene.pos$Start, gene.pos$End)
gene.pos$End.dir <- ifelse(gene.pos$Strand=='+', gene.pos$End, gene.pos$Start)
gene.pos$midpoint <- as.numeric(gene.pos$Start) + ((as.numeric(gene.pos$End) - as.numeric(gene.pos$Start) )/2)
genetext <- gene.pos$Gene

library(gggenes)
p.genepos <- ggplot(gene.pos,aes(xmin = Start.dir, xmax = End.dir, y = "molecule", fill = Strand)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm")) + 
  scale_x_continuous(labels = scales::comma,expand = c(0, 0), breaks=seq(0, genomelength, round(genomelength / 15,1 -nchar(ceiling(genomelength / 15))))) +
  coord_cartesian(xlim=c(0,genomelength)) +
  theme_bw() + theme(legend.position="none") +
  theme( axis.text.y = element_blank(), axis.ticks.y= element_blank()) + labs(y="Genes")
p.genepos <- p.genepos + geom_text(data=gene.pos, aes(x=midpoint, y="Genes", label=Gene),angle = -90, size=2)


y.theme.strip <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y= element_blank())
x.theme.strip <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x= element_blank())

g.genepos <- ggplotGrob(p.genepos + x.theme.strip)
g.SNPdensity <- ggplotGrob(p.SNPdensity)
g.genepos$widths <- g.SNPdensity$widths

grid.arrange(g.SNPdensity,g.genepos,ncol=1,heights=c(12,1))

 
Cairo(file=paste(genefile,".SNP-density-facetted.png", sep="") , width = 1000, height = 1000,type="png",dpi=600, units = "pt")
#p.SNPdensity
grid.arrange(g.SNPdensity,g.genepos,ncol=1,heights=c(12,1))
dev.off()



###




  
# Generate png file of plot (in same location as original fasta)
Cairo(file=paste(genefile,".SNP-plot.png", sep="") , width = 900, height = 700,type="png",dpi=600, units = "pt")
p2
dev.off()
# Generate a tab-delimited file of variant presence/absence 
write.table(gene.aln.binary2, file=paste(genefile,".variant-matrix.tdl", sep=""),quote=F, sep="\t", col.names=T )



####
# Plot SNP heatmap against a tree
# N.B. Tip labels must be identical to sample samples in sequence alignment
library(ggtree)


#First bring in the tree

#treefile <- "/Users/mb29/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Gubbins/Szpara2014-v4/raxml-bootstrap/RAxML_bipartitionsBranchLabels.HSV1_Corededup+Szpara2014_v3.tre"
treefile <- "/Users/mb29/Treponema/Lukehart_UoW/MiniProjects1_intrahost/Analysis/Phylo/iqtree_WGS/SLK3a_bwa.aln.treefile"
mytree <- read.tree(treefile)  
mytree1 <- ggtree(mytree) 
mytree1 <- ggtree(phytools::midpoint.root(mytree))
mytree1

# Reorder the dataframe/heatmap
gene.aln.windowed.decast <-dcast(gene.aln.windowed.melt, Sequence~Position, value.var="Variant")
rownames(gene.aln.windowed.decast) <- gene.aln.windowed.decast[,1]
gene.aln.windowed.decast <- gene.aln.windowed.decast[,c(2:ncol(gene.aln.windowed.decast))]

# Plot both together (might want to tweak the colours and options)
mytree1 %>% gheatmap(gene.aln.windowed.decast,width=3,font=3,color=NULL) +
  #scale_color_manual(breaks=c("Reference", "Variant"),values=c("white", "#F8766D")) +
  scale_fill_manual(breaks=c("Reference", "Variant"),values=c("white", "#F8766D")) + geom_tiplab()



##########
# PCA on binary matrix
gene.aln.binary.mat <- gene.aln.binary[,c(1:(ncol(gene.aln.binary)-1))]
row.names(gene.aln.binary.mat) <- gene.aln.binary$position

gene.aln.binary.pca <- prcomp(t(data.matrix(gene.aln.binary.mat, rownames.force = NA)),center = T)
gene.aln.binary.pca <- as.data.frame(gene.aln.binary.pca$x)
gene.aln.binary.pca $Sample <- rownames(gene.aln.binary.pca )
gene.aln.binary.pca$SampleSet <- ifelse(grepl("SWAB",gene.aln.binary.pca$Sample),"SWAB","CSF")
gene.aln.binary.pca$SampleSet <- sapply(1:nrow(gene.aln.binary.pca), function(x) ifelse(grepl("NC001806.2_HSV1_s_17",gene.aln.binary.pca$Sample[x]),"Reference",gene.aln.binary.pca$SampleSet[x]))

gene.aln.binary.pca$Sample <- gsub("\\_HSV1","",gsub("HSV1-","",gene.aln.binary.pca$Sample))

p.pca <- ggplot(gene.aln.binary.pca,aes(PC1,PC2, color=factor(SampleSet),label=Sample)) + 
  #geom_point(alpha=0.5) +
  geom_text(size=3.5) + 
  theme_bw() +
  scale_colour_manual(values = c("#009900","black","darkorchid3")) +
  theme(legend.position="top") + labs(color="Sample Set") 
p.pca

Cairo(file=paste(genefile,".SNPs-PCA.png", sep="") , width = 600, height = 500,type="png",dpi=600, units = "pt")
p.pca
dev.off()


# plot tSNE
library(Rtsne)
gene.aln.binary.tsne <- Rtsne(t(data.matrix(gene.aln.binary.mat, rownames.force = NA)), dims=2, perplexity=3, verbose=TRUE, max_iter = 500)
gene.aln.binary.tsne <- as.data.frame(gene.aln.binary.tsne$Y)
gene.aln.binary.tsne$Sample <- colnames(gene.aln.binary.mat)
gene.aln.binary.tsne$SampleSet <- ifelse(grepl("SWAB",gene.aln.binary.tsne$Sample),"SWAB","CSF")
gene.aln.binary.tsne$SampleSet <- sapply(1:nrow(gene.aln.binary.tsne), function(x) ifelse(grepl("NC001806.2_HSV1_s_17",gene.aln.binary.tsne$Sample[x]),"Reference",gene.aln.binary.tsne$SampleSet[x]))

p.tsne <- ggplot(gene.aln.binary.tsne,aes(V1,V2, color=factor(SampleSet),label=Sample)) + 
  #geom_point(alpha=0.5) +
  geom_text(size=3.5) + 
  theme_bw() +
  scale_colour_manual(values = c("#009900","black","darkorchid3")) +
  theme(legend.position="top") + labs(color="Sample Set") 
p.pca
