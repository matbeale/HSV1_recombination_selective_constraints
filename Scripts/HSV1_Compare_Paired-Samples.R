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
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
#biocLite("GenomicFeatures")
#biocLite("BSgenome")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("BSgenome")
#BiocManager::install("limma")

library("VariantAnnotation")
library("GenomicFeatures")
library("BSgenome")

linked.vcf.dir <- "/Users/mb29/hsv1/freebayes-vcfs/Linked_Pairs/"


Pair1 <- c("HSV1-nCSF1","HSV1-nCSF3")
Pair2 <- c("HSV1-nCSF10","HSV1-nCSF9")
Pair3 <- c("HSV1-CSF4","HSV1-CSF7","HSV1-CSF8")
Pair4 <- c("HSV1-CSF5","HSV1-CSF6")


genomelength <- 152222
myvar.collating <- data.frame(c(1:genomelength),stringsAsFactors = F)
colnames(myvar.collating) <- "position"

for (currentseq in Pair3) {

  #vcf <- readVcf(paste("C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/MinorVars/Linked_Pairs/",currentseq,"-freebayes.5-reads.vcf",sep=""),"NC_001806.2_HSV1_s17")
  vcf <- readVcf(paste(linked.vcf.dir,currentseq,"-freebayes.5-reads.vcf",sep=""),"NC_001806.2_HSV1_s17")
  
  myvcf <- data.frame(info(vcf))
  myvcf$position <- as.numeric(gsub("\\_.+$","",gsub("^.+\\:","", rownames(myvcf),perl=T), perl=T))
  myvcf <- myvcf[,c("position","NS","DP","AN","RO","AO","SRF","SRR","SAF","SAR")]
  colnames(myvcf)[c(3:6)] <- c("Depth","Alleles","RefCount","AltCount")
  myvcf$keep <- sapply(1:nrow(myvcf), function(x) ifelse((grepl(",",myvcf$AltCount[x]) ),"drop","keep"))
  myvcf$keep <- sapply(1:nrow(myvcf), function(x) ifelse((grepl(",",myvcf$RefCount[x]) ),"drop",myvcf$keep[x]))
  myvcf <- subset(myvcf,keep=="keep", select=c(1:ncol(myvcf)-1))
  myvcf$Depth2 <- sapply(1:nrow(myvcf), function(x) sum(as.numeric(myvcf$RefCount[x]), as.numeric(myvcf$AltCount[x]))) # depth provided in VCF is missing some reads, so recalculate this
  myvcf$refper <- round(((as.numeric(myvcf$RefCount)/as.numeric(myvcf$Depth2))*100),1)
  myvcf$altper <- round(((as.numeric(myvcf$AltCount)/as.numeric(myvcf$Depth2))*100),1)
  # Produce filter to keep if >50% alt variant
  myvcf.minors <- myvcf # keep a backup for pulling out minor vars
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(myvcf$altper[x] >=50, "alt","ref" ))
  myvcf <- subset(myvcf, vartype=="alt")
  # Modify filter to remove minor vars with low coverage
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$AltCount[x]) <5, "drop",myvcf$vartype[x]))
  myvcf <- subset(myvcf, vartype=="alt")
  # Add in filter to remove variants not covered by at least one read on each strand
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$SAF[x]) <1, "drop",myvcf$vartype[x]))
  myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$SAR[x]) <1, "drop",myvcf$vartype[x]))
  myvcf <- subset(myvcf, vartype=="alt")
  
  # Go back to earlier saved version and produce filter to remove "major variant sites"
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(myvcf.minors$refper[x] >= 98 || myvcf.minors$altper[x] >=98, "major","minor" ))
  myvcf.minors <- subset(myvcf.minors, vartype=="minor")
  # Modify filter to remove minor vars with low coverage
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$RefCount[x]) <5, "drop",myvcf.minors$vartype[x]))
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$AltCount[x]) <5, "drop",myvcf.minors$vartype[x]))
  myvcf.minors <- subset(myvcf.minors, vartype=="minor")
  myvcf.minors$majorfreq <- sapply(1:nrow(myvcf.minors), function (x) max(c(as.numeric(myvcf.minors$refper[x]) ,as.numeric(myvcf.minors$altper[x]))))
  myvcf.minors$minorfreq <- sapply(1:nrow(myvcf.minors), function (x) min(c(as.numeric(myvcf.minors$refper[x]) ,as.numeric(myvcf.minors$altper[x]))))
    # Add in filter to remove minor variants not covered by at least one read on each strand
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SRF[x]) <1, "drop",myvcf.minors$vartype[x]))
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SRR[x]) <1, "drop",myvcf.minors$vartype[x]))
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SAF[x]) <1, "drop",myvcf.minors$vartype[x]))
  myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SAR[x]) <1, "drop",myvcf.minors$vartype[x]))
  myvcf.minors <- subset(myvcf.minors, vartype=="minor")
  
  # Add test of strand bias, and remove variants that that fail test 
  #strandbias.cuttoff <- 0.1 
  #myvcf.minors$RefStrBias <- sapply(1:nrow(myvcf.minors), function(x) round(as.numeric(as.numeric(myvcf.minors$SRF[x])/sum(as.numeric(myvcf.minors$SRF[x]),as.numeric(myvcf.minors$SRR[x]))),2))
  #myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$RefStrBias[x]) < strandbias.cuttoff || as.numeric(myvcf.minors$RefStrBias[x]) >1-strandbias.cuttoff, "drop",myvcf.minors$vartype[x]))
  #myvcf.minors$AltStrBias <- sapply(1:nrow(myvcf.minors), function(x) round(as.numeric(as.numeric(myvcf.minors$SAF[x])/sum(as.numeric(myvcf.minors$SAF[x]),as.numeric(myvcf.minors$SAR[x]))),2))
  #myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$AltStrBias[x]) < strandbias.cuttoff || as.numeric(myvcf.minors$AltStrBias[x]) >>1-strandbias.cuttoff, "drop",myvcf.minors$vartype[x]))
  #myvcf.minors <- subset(myvcf.minors, vartype=="minor", select=-c(RefStrBias,AltStrBias))
  
  # Now bring Major and Minor Vars back together
  myvar.called <- data.frame(myvcf[,c("position","altper")],stringsAsFactors = F)
  colnames(myvar.called)[ncol(myvar.called)] <- currentseq
  myvar.collating <- merge(myvar.collating,myvar.called, by.x="position",by.y="position",all.x=T)
  myvar.collating[ncol(myvar.collating)] <- sapply(1:nrow(myvar.collating), function(x) ifelse(is.na(myvar.collating[x,ncol(myvar.collating)]),0,1))
  
  myvcf.minors.called <- data.frame(myvcf.minors[,c("position","altper")],stringsAsFactors = F)
  colnames(myvcf.minors.called)[ncol(myvcf.minors.called)] <- paste(currentseq,".altper", sep="")
  myvar.collating <- merge(myvar.collating,myvcf.minors.called, by.x="position",by.y="position",all.x=T)
}

# strip out positions without any variants in dataset
myvar.collating[is.na(myvar.collating)]<- 0
myvar.collating$keep <- sapply(1:nrow(myvar.collating), function(x) sum(myvar.collating[x,2:ncol(myvar.collating)]))
myvar.collating <- subset(myvar.collating, keep>0, select=-keep)
nrow(myvar.collating)

# Identify any sites where the consensus changes
myvar.cons.change <- myvar.collating
myvar.cons.change$keep <- sapply(1:nrow(myvar.cons.change), function (x) sum(myvar.cons.change[x,grep("altper",colnames(myvar.cons.change), invert=T)[-1]]))
myvar.cons.change <- subset(myvar.cons.change, keep==1, select=-keep)
#myvar.cons.change <- subset(myvar.cons.change, keep!=3, select=-keep)
#myvar.cons.change <- subset(myvar.cons.change, keep!=0, select=-keep)
nrow(myvar.cons.change)

# Keep only sites with minor variants
myvar.minor.change <- myvar.collating
myvar.minor.change$keep <- sapply(1:nrow(myvar.minor.change), function (x) sum(myvar.minor.change[x,grep("altper",colnames(myvar.minor.change))]))
myvar.minor.change <- subset(myvar.minor.change, keep!=0, select=-keep)

# Now try and plot
myvar.minor.melt <- myvar.minor.change[,c(1,grep("altper",colnames(myvar.minor.change)))]
colnames(myvar.minor.melt) <- gsub(".altper","",colnames(myvar.minor.melt))
myvar.minor.melt <- melt(myvar.minor.melt, id.vars="position")
myvar.minor.melt$position <- as.factor(myvar.minor.melt$position)
# rename nCSF as SWAB
myvar.minor.melt$variable <- gsub("nCSF","SWAB",myvar.minor.melt$variable)

p.minor <- ggplot(myvar.minor.melt, aes(position,value, group=variable,  fill=variable))
p.minor <- p.minor + geom_bar(stat="identity") + theme_bw() +
  facet_grid(variable ~ ., scales="free", space="free", drop=T, margins=F) +
  labs(x="Variant position",y="Alt Variant Frequency (%)") + coord_cartesian(ylim=c(2,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5), legend.position="none", 
        strip.background = element_rect(fill="white"))


Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair3","_MinorVarSites.svg", sep=""), width = 700, height = 300, type="svg", units = "pt")
#Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair1","_MinorVarSites.png", sep=""), width = 700, height = 300, type="png",dpi=600, units = "pt")
p.minor
dev.off()



library("VennDiagram")
#library("venneuler")
library("limma")
library("extrafont")
#font_import()

# Pull out minor var sites
myvar.minor.venn <- myvar.minor.change[,c(1,grep("altper",colnames(myvar.minor.change)))]
myvar.minor.venn[myvar.minor.venn > 0] <- 1
myvar.minor.venn <- data.frame(myvar.minor.venn[2:ncol(myvar.minor.venn)], stringsAsFactors = F)
colnames(myvar.minor.venn) <- gsub("\\.","-",(gsub(".altper","",colnames(myvar.minor.venn))))
myvar.minor.vcounts <- vennCounts(myvar.minor.venn)
# Relabel nCSF as SWAB
colnames(myvar.minor.vcounts) <- gsub("nCSF","SWAB",colnames(myvar.minor.vcounts))

# 2-way Venn
grid.newpage()
myvenn <- draw.pairwise.venn(area2 = sum(myvar.minor.vcounts[2,3],myvar.minor.vcounts[4,3]), 
                             area1 = sum(myvar.minor.vcounts[3,3],myvar.minor.vcounts[4,3]), 
                             cross.area = myvar.minor.vcounts[4,3], 
                             category = gsub("HSV1-","",colnames(myvar.minor.vcounts)[1:2]), 
                             fill = c("coral","cyan"),alpha = c(0.75,0.4), cex = 2,  cat.cex = 1.5,
                             rotation.degree=-90, cat.dist = 0.025,margin = 0.05,cat.pos=c(0,180),
                             #cex.prop=2,#inverted=T,
                             fontfamily = "Arial", main.fontfamily="Arial",cat.fontfamily="Arial",
                             scaled = T, lty = rep("blank", 2),main="Minor Variant Sites between Pair2", ind=F)
grid.draw(myvenn)


# 3-way Venn
grid.newpage()
myvenn <- draw.triple.venn(area1= sum(myvar.minor.vcounts[5,4],myvar.minor.vcounts[6,4],myvar.minor.vcounts[7,4],myvar.minor.vcounts[8,4]) ,
                           area2= sum(myvar.minor.vcounts[3,4],myvar.minor.vcounts[4,4],myvar.minor.vcounts[7,4],myvar.minor.vcounts[8,4]), 
                           area3= sum(myvar.minor.vcounts[2,4],myvar.minor.vcounts[4,4],myvar.minor.vcounts[6,4],myvar.minor.vcounts[8,4]) ,
                           n12= sum(myvar.minor.vcounts[7,4],myvar.minor.vcounts[8,4]), n23= sum(myvar.minor.vcounts[4,4],myvar.minor.vcounts[8,4]),n13=sum(myvar.minor.vcounts[6,4],myvar.minor.vcounts[8,4]),n123=myvar.minor.vcounts[8,4],
                           category = gsub("HSV1-","",colnames(myvar.minor.vcounts)[1:3]),
                           fill = c("coral", "green1","cyan"), alpha = c(0.6,0.3, 0.7), cex = 2,  cat.cex = 1.5,
                           rotation.degree=90, cat.dist = 0.075,margin = 0.05,cat.pos=c(180,0,270),#cat.pos=c(300,60,200),
                           fontfamily = "Arial", main.fontfamily="Arial",cat.fontfamily="Arial",
                           scaled = T, lty = rep("blank", 3), ind=F)
grid.draw(myvenn)                           
                            

Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair3","_MinorVarSites_Venn.svg", sep=""), width = 430, height = 430, type="svg", units = "pt")
#Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair2","_MinorVarSites_Venn.png", sep=""), width = 430, height = 430,type="png", units = "pt",dpi=900)
grid.draw(myvenn)
#grid.newpage()
#draw.pairwise.venn(area1 = sum(myvar.minor.vcounts[2,3],myvar.minor.vcounts[4,3]), area2 = sum(myvar.minor.vcounts[3,3],myvar.minor.vcounts[4,3]), cross.area = myvar.minor.vcounts[4,3], category = colnames(myvar.minor.vcounts)[1:2], fill = c("coral","cyan"),alpha = c(0.75,0.4),scaled = T, lty = rep("blank", 2))
dev.off()




########################
#
#
#
################ Do more detailed evaluation of Pair 1 - look at relative coverage, etc.

#highcov.names <- c("HSV1-nCSF11","HSV1-nCSF10","HSV1-nCSF7","HSV1-nCSF14","HSV1-nCSF15","HSV1-nCSF5","HSV1-CSF1","HSV1-nCSF1","HSV1-CSF4","HSV1-CSF2","HSV1-CSF5","HSV1-nCSF6","HSV1-CSF3","HSV1-nCSF13")

highcov.names <- Pair1
coverage.threshold <- 1000
sample.threshold <- length(highcov.names)

minorvars.cov <- data.frame(myvar.minor.change$position, stringsAsFactors = F)
colnames(minorvars.cov) <- "position"

for (currentseq in highcov.names) {  
  covpath <- linked.vcf.dir
  seqname <- paste(currentseq,".depth", sep="")
  coverage.called <- read.table(paste(covpath,seqname,sep=""), header=T, stringsAsFactors = F)
  colnames(coverage.called) <- c("ref","position",paste(currentseq,".cov",sep=""))
  coverage.called <- coverage.called[,-1]
  # collate coverage at each minor site
  minorvars.cov <- merge(minorvars.cov, coverage.called, by="position")
  # evaluate if site meets minimum threshold (250x)
  minorvars.cov[ncol(minorvars.cov)] <- sapply(1:nrow(minorvars.cov), function(x) ifelse(minorvars.cov[x,ncol(minorvars.cov)] <coverage.threshold,0,1))
}
minorvars.cov$passed <- sapply(1:nrow(minorvars.cov), function(x) sum(minorvars.cov[x,c(2:ncol(minorvars.cov))]))

# Only use highcov
highcov.minorvars <- myvar.minor.change[,c(1,3,5)] #myvar.minor.change[,colnames(myvar.minor.change)%in%highcov.names]
#highcov.minorvars[highcov.minorvars>=2]<- 1 # Convert frequency data into binary value
#highcov.minorvars[highcov.minorvars<2]<- 0 # Ensure any low values are made 0 (should already be done)
highcov.minorvars$position <- rownames(highcov.minorvars)
minorvars.cov.passed <- minorvars.cov[,c(1,4)]
highcov.minorvars <- merge(highcov.minorvars, minorvars.cov.passed, by="position")
nrow(subset(highcov.minorvars, passed==sample.threshold, select=-passed)) # get a count
highcov.minorvars <- subset(highcov.minorvars, passed==sample.threshold, select=-passed) # subset to get highcov sites

# Now try and plot
#highcov.minor.melt <- highcov.minorvars[,c(1,grep("altper",colnames(highcov.minorvars)))]
highcov.minor.melt <- highcov.minorvars
colnames(highcov.minor.melt) <- gsub(".altper","",colnames(highcov.minor.melt))
highcov.minor.melt <- melt(highcov.minor.melt, id.vars="position")
highcov.minor.melt$position <- as.factor(highcov.minor.melt$position)

p.high.minor <- ggplot(highcov.minor.melt, aes(position,value, group=variable,  fill=variable))
p.high.minor <- p.high.minor + geom_bar(stat="identity") + theme_bw() +
  facet_grid(variable ~ ., scales="free", space="free", drop=T, margins=F) +
  labs(x="Variant position",y="Alt Variant Frequency (%)") + coord_cartesian(ylim=c(2,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5), legend.position="none")
p.high.minor


Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair1","_MinorVarSites_pass+",coverage.threshold,"x-cov.svg", sep=""), width = 800, height = 400, type="svg", units = "pt")
#Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair1","_MinorVarSites_pass+",coverage.threshold,"x-cov.png", sep=""), width = 800, height = 400, type="png",dpi=600, units = "pt")
p.high.minor
dev.off()

highcov.minorvars.venn <- highcov.minorvars#myvar.minor.change[,c(1,grep("altper",colnames(myvar.minor.change)))]
highcov.minorvars.venn[highcov.minorvars.venn > 0] <- 1
highcov.minorvars.venn <- data.frame(highcov.minorvars.venn[2:ncol(highcov.minorvars.venn)], stringsAsFactors = F)
colnames(highcov.minorvars.venn) <- gsub("\\.","-",(gsub(".altper","",colnames(highcov.minorvars.venn))))
highcov.minorvars.vcounts <- vennCounts(highcov.minorvars.venn)

# 2-way Venn
grid.newpage()
myvenn <- draw.pairwise.venn(area2 = sum(highcov.minorvars.vcounts[2,3],highcov.minorvars.vcounts[4,3]), 
                             area1 = sum(highcov.minorvars.vcounts[3,3],highcov.minorvars.vcounts[4,3]), 
                             cross.area = highcov.minorvars.vcounts[4,3], 
                             category = colnames(highcov.minorvars.vcounts)[1:2], 
                             fill = c("coral","cyan"),alpha = c(0.75,0.4),
                             scaled = T, lty = rep("blank", 2),main="Minor Variant Sites between Pair1", ind=F,fontfamily="Arial")
grid.draw(myvenn)


Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair1","_MinorVarSites_pass+",coverage.threshold,"x-cov_Venn.svg", sep=""), width = 430, height = 430, type="svg", units = "pt")
#Cairo(file= paste(linked.vcf.dir,"HSV1_","Pair1","_MinorVarSites_pass+",coverage.threshold,"x-cov_Venn.png", sep=""), width = 430, height = 430,type="png", units = "pt",dpi=900)
grid.draw(myvenn)
dev.off()





###################
#
#
# Look at VCFs manually
#Pair1 <- c("HSV1-nCSF1","HSV1-nCSF3")
Pair1.1 <- "HSV1-nCSF1"
Pair1.2 <- "HSV1-nCSF3"

vcf <- readVcf(paste(linked.vcf.dir,Pair1.1,"-freebayes.5-reads.vcf",sep=""),"NC_001806.2_HSV1_s17")
#vcf <- readVcf(paste("C:/Bioinformatics/HSV1_genomes/annotation/snippy_0.5-cutoff/MinorVars/Linked_Pairs/",Pair1.2,"-freebayes.5-reads.vcf",sep=""),"NC_001806.2_HSV1_s17")



myvcf <- data.frame(info(vcf))
myvcf$position <- as.numeric(gsub("\\_.+$","",gsub("^.+\\:","", rownames(myvcf),perl=T), perl=T))
myvcf <- myvcf[,c("position","NS","DP","AN","RO","AO","SRF","SRR","SAF","SAR")]
colnames(myvcf)[c(3:6)] <- c("Depth","Alleles","RefCount","AltCount")
myvcf$keep <- sapply(1:nrow(myvcf), function(x) ifelse((grepl(",",myvcf$AltCount[x]) ),"drop","keep"))
myvcf$keep <- sapply(1:nrow(myvcf), function(x) ifelse((grepl(",",myvcf$RefCount[x]) ),"drop",myvcf$keep[x]))
myvcf <- subset(myvcf,keep=="keep", select=c(1:ncol(myvcf)-1))
myvcf$Depth2 <- sapply(1:nrow(myvcf), function(x) sum(as.numeric(myvcf$RefCount[x]), as.numeric(myvcf$AltCount[x]))) # depth provided in VCF is missing some reads, so recalculate this
myvcf$refper <- round(((as.numeric(myvcf$RefCount)/as.numeric(myvcf$Depth2))*100),1)
myvcf$altper <- round(((as.numeric(myvcf$AltCount)/as.numeric(myvcf$Depth2))*100),1)
# Produce filter to keep if >50% alt variant
myvcf.minors <- myvcf # keep a backup for pulling out minor vars
myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(myvcf$altper[x] >=50, "alt","ref" ))
myvcf <- subset(myvcf, vartype=="alt")
# Modify filter to remove minor vars with low coverage
myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$AltCount[x]) <5, "drop",myvcf$vartype[x]))
myvcf <- subset(myvcf, vartype=="alt")
# Add in filter to remove variants not covered by at least one read on each strand
myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$SAF[x]) <1, "drop",myvcf$vartype[x]))
myvcf$vartype <- sapply(1:nrow(myvcf), function(x) ifelse(as.numeric(myvcf$SAR[x]) <1, "drop",myvcf$vartype[x]))
myvcf <- subset(myvcf, vartype=="alt")

# Go back to earlier saved version and produce filter to remove "major variant sites"
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(myvcf.minors$refper[x] >= 98 || myvcf.minors$altper[x] >=98, "major","minor" ))
myvcf.minors <- subset(myvcf.minors, vartype=="minor")
# Modify filter to remove minor vars with low coverage
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$RefCount[x]) <5, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$AltCount[x]) <5, "drop",myvcf.minors$vartype[x]))
myvcf.minors <- subset(myvcf.minors, vartype=="minor")
myvcf.minors$majorfreq <- sapply(1:nrow(myvcf.minors), function (x) max(c(as.numeric(myvcf.minors$refper[x]) ,as.numeric(myvcf.minors$altper[x]))))
myvcf.minors$minorfreq <- sapply(1:nrow(myvcf.minors), function (x) min(c(as.numeric(myvcf.minors$refper[x]) ,as.numeric(myvcf.minors$altper[x]))))
# Add in filter to remove minor variants not covered by at least one read on each strand
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SRF[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SRR[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SAF[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$SAR[x]) <1, "drop",myvcf.minors$vartype[x]))
myvcf.minors <- subset(myvcf.minors, vartype=="minor")

# Add test of strand bias, and remove variants that that fail test
strandbias.cuttoff <- 0.1
myvcf.minors$RefStrBias <- sapply(1:nrow(myvcf.minors), function(x) round(as.numeric(as.numeric(myvcf.minors$SRF[x])/sum(as.numeric(myvcf.minors$SRF[x]),as.numeric(myvcf.minors$SRR[x]))),2))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$RefStrBias[x]) < strandbias.cuttoff || as.numeric(myvcf.minors$RefStrBias[x]) >1-strandbias.cuttoff, "drop",myvcf.minors$vartype[x]))
myvcf.minors$AltStrBias <- sapply(1:nrow(myvcf.minors), function(x) round(as.numeric(as.numeric(myvcf.minors$SAF[x])/sum(as.numeric(myvcf.minors$SAF[x]),as.numeric(myvcf.minors$SAR[x]))),2))
myvcf.minors$vartype <- sapply(1:nrow(myvcf.minors), function(x) ifelse(as.numeric(myvcf.minors$AltStrBias[x]) < strandbias.cuttoff || as.numeric(myvcf.minors$AltStrBias[x]) >1-strandbias.cuttoff, "drop",myvcf.minors$vartype[x]))
#myvcf.minors <- subset(myvcf.minors, vartype=="minor", select=-c(RefStrBias,AltStrBias))
#myvcf.minors <- subset(myvcf.minors, vartype=="minor")









