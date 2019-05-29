
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

####

seqpath <- "C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/VCFs/"

myvar.collating <- read.table(paste(seqpath,"HSV1-All-MajorVar-sites_merged.freebayes.txt", sep=""), stringsAsFactors = F, header=T)
colnames(myvar.collating) <- gsub("\\.","\\-",colnames(myvar.collating))

# strip out positions without any variants in dataset
myvar.collating$keep <- sapply(1:nrow(myvar.collating), function(x) sum(myvar.collating[x,2:ncol(myvar.collating)]))
myvar.collating <- subset(myvar.collating, keep>0, select=-keep)
nrow(myvar.collating)


wgsaln.file <- "C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/Core_Snps/HSV1_Core-dedup.full.fas"
wgsaln <- seqinr::read.fasta(file=wgsaln.file, seqtype="DNA", as.string=T,set.attributes=F)
seqnames <- names(wgsaln)
wgsaln <- data.frame(t(as.matrix.alignment(as.alignment(nb=19,nam=NULL,seq=wgsaln,com=NULL))), stringsAsFactors = F)
colnames(wgsaln) <- seqnames
wgsaln$position <- rownames(wgsaln)

wgsaln.sites <- wgsaln[rownames(wgsaln)%in%myvar.collating$position,]
            




###############
## Try applying a Fisher's exact test to variants (GWAS) between populations
#####
snpmatrix3 <- wgsaln.sites
#snpmatrix3 <- merge(snpmatrix, genes.snp.no, by="snpmatrix")

# Determine number of alleles at each site
snpmatrix3$alleles <- sapply(1:nrow(snpmatrix3), function(x) length(unique(as.character(snpmatrix3[x,c(2:19)]))))
# filter for biallelic sites only (initially)
nrow(subset(snpmatrix3, alleles == 3)) # 22 triallelic sites (185 including 'missing')
nrow(subset(snpmatrix3, alleles == 2)) # 2562 biallelic sites in this 'all sites' dataset (not from gene based)
snpmatrix3.biallelic <- subset(snpmatrix3, alleles == 2, select=-c(20:21))
snps.alleles <- data.frame(rownames(snpmatrix3.biallelic), stringsAsFactors = F)
colnames(snps.alleles) <- "position"
snps.alleles$var1 <- sapply(1:nrow(snpmatrix3.biallelic), function(x) unique(as.character(snpmatrix3.biallelic[x,c(2:19)]))[1])
snps.alleles$var2 <- sapply(1:nrow(snpmatrix3.biallelic), function(x) unique(as.character(snpmatrix3.biallelic[x,c(2:19)]))[2])
# subset for CSF and determine frequency of  each allele
snpsmatrix.transposed <- data.frame(t(snpmatrix3.biallelic),stringsAsFactors = F)
snpsmatrix.transposed$seqs <- rownames(snpsmatrix.transposed)
snpsmatrix.transposed <- merge(snpsmatrix.transposed, Sampletype, by.x="seqs", by.y="Sample")
snpsmatrix.transposed <- subset(snpsmatrix.transposed, Source=="CSF")#, select=c(3:ncol(snps.transposed)-1))
snps.alleles$CSFvar1 <- sapply(1:nrow(snps.alleles), function(x) length(snpsmatrix.transposed[snpsmatrix.transposed[x+1]==snps.alleles$var1[x],x+1]))
snps.alleles$CSFvar2 <- sapply(1:nrow(snps.alleles), function(x) length(snpsmatrix.transposed[snpsmatrix.transposed[x+1]==snps.alleles$var2[x],x+1]))
# subset for non-CSF and determine frequency of  each allele
snpsmatrix.transposed <- data.frame(t(snpmatrix3.biallelic),stringsAsFactors = F)
snpsmatrix.transposed$seqs <- rownames(snpsmatrix.transposed)
snpsmatrix.transposed <- merge(snpsmatrix.transposed, Sampletype, by.x="seqs", by.y="Sample")
snpsmatrix.transposed <- subset(snpsmatrix.transposed, Source=="nonCSF")#, select=c(3:ncol(snps.transposed)-1))
snps.alleles$nCSFvar1 <- sapply(1:nrow(snps.alleles), function(x) length(snpsmatrix.transposed[snpsmatrix.transposed[x+1]==snps.alleles$var1[x],x+1]))
snps.alleles$nCSFvar2 <- sapply(1:nrow(snps.alleles), function(x) length(snpsmatrix.transposed[snpsmatrix.transposed[x+1]==snps.alleles$var2[x],x+1]))

# bring in "consequence data"
#snps.alleles$gene <- gsub("\\_.+$","",snps.alleles$snpmatrix3.biallelic.snpmatrix)
#snps.alleles <- merge(snps.alleles, genes.snp.no, by.x="snpmatrix3.biallelic.snpmatrix" ,by.y="snpmatrix")

################
# Filter to only test genes of interest (from selection analysis)
#snps.alleles$keep <- ifelse((snps.alleles$gene == "UL18" | snps.alleles$gene == "UL33" | snps.alleles$gene == "UL51"), "keep","drop")
#snps.alleles <- subset(snps.alleles, keep == "keep", select=-11)

# Filter to only test non-syn snps
#snps.alleles <- subset(snps.alleles, syn == "non-syn")
################

# create 2x2 contingency table for each loci and test using Fisher's exact test
snps.alleles$fisher <- unlist(sapply(1:nrow(snps.alleles), function(x) fisher.test(matrix(c(snps.alleles$CSFvar1[x], snps.alleles$CSFvar2[x], snps.alleles$nCSFvar1[x],snps.alleles$nCSFvar2[x]),nrow=2))[1]))
# correct for multiple testing
snps.alleles$bonferroni <- min(p.adjust(snps.alleles$fisher, method="bonferroni"))
snps.alleles$fdr <- min(p.adjust(snps.alleles$fisher, method="fdr"))
#### And we find that no loci are significant after adjustment for multiple testing ####
##################
##################







