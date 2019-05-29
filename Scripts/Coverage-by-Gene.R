#!/usr/bin/env Rscript

#Usage : HBV-WGS_2_Resistance_variants_CL.R \
#  <Consensus.fasta>   # 1
#  <annotationsfile>   # 2
#  <varscan2file.txt>  # 3
#  <varscan20.txt>     # 4 
#  <variantlist>       # 5
#  <refseqfile>        # 6
#  <bamfile>           # 7
###  <coveragefile>      # 8

coveragefile <- "C:/Users/Mat/Dropbox/UCL_Pathseek/MyScripts/Local_R/varscan/HBV-2116.perbase.coverage"
positionfile <- "C:/Users/Mat/Dropbox/UCL_Pathseek/MyScripts/Local_R/varscan/HBV-2116.annotations"


#args <- commandArgs(trailingOnly = TRUE) # capture command line arguments
# Put command line arguments into variables
#consensusfile <- args[1]
#positionfile <- args[2]
#varscan2file <- args[3]
#varscan20file <- args[4]
#mutantfile <- args[5]
#reffile <- args[6]
#bamFile <- args[7]
#coveragefile <- args[8]

coverage = read.table(coveragefile,header=F)  #Read in file
positions <- read.table(positionfile, header=F)
seqname <- sub(".annotations","", positionfile, ignore.case=TRUE, perl=TRUE)
sequence <- sub("^.+HBV-", "HBV-",seqname, perl=T)

# Whole Genome stats
WGS.100x <- (sum(coverage$V3 >=100)/nrow(coverage))*100
WGS.2x <- (sum(coverage$V3 >1)/nrow(coverage))*100
WGS.mean <- round(mean(coverage$V3),1)
WGS.min  <- min(coverage$V3)

# Polymerase
Polcoverage <- coverage[c(positions[7,2]:positions[7,3],positions[1,2]:positions[1,3]),]
Pol.100x <- (sum(Polcoverage$V3 >=100)/nrow(Polcoverage))*100
Pol.2x <- (sum(Polcoverage$V3 >1)/nrow(Polcoverage))*100
Pol.mean <- round(mean(Polcoverage$V3),1)
Pol.min  <- min(Polcoverage$V3)

# RT
RTcoverage <- coverage[c(positions[2,2]:positions[2,3]),]
RT.100x <- (sum(RTcoverage$V3 >=100)/nrow(RTcoverage))*100
RT.2x <- (sum(RTcoverage$V3 >1)/nrow(RTcoverage))*100
RT.mean <- round(mean(RTcoverage$V3),1)
RT.min  <- min(RTcoverage$V3)

# Surface
Surfacecoverage <- coverage[c(positions[3,2]:positions[3,3]),]
Surface.100x <- (sum(Surfacecoverage$V3 >=100)/nrow(Surfacecoverage))*100
Surface.2x <- (sum(Surfacecoverage$V3 >1)/nrow(Surfacecoverage))*100
Surface.mean <- round(mean(Surfacecoverage$V3),1)
Surface.min  <- min(Surfacecoverage$V3)

# X
Xcoverage <- coverage[c(positions[4,2]:positions[4,3]),]
X.100x <- (sum(Xcoverage$V3 >=100)/nrow(Xcoverage))*100
X.2x <- (sum(Xcoverage$V3 >1)/nrow(Xcoverage))*100
X.mean <- round(mean(Xcoverage$V3),1)
X.min  <- min(Xcoverage$V3)

# PreC
PreCcoverage <- coverage[c(positions[5,2]:positions[5,3]),]
PreC.100x <- (sum(PreCcoverage$V3 >=100)/nrow(PreCcoverage))*100
PreC.2x <- (sum(PreCcoverage$V3 >1)/nrow(PreCcoverage))*100
PreC.mean <- round(mean(PreCcoverage$V3),1)
PreC.min  <- min(PreCcoverage$V3)

# Core
Corecoverage <- coverage[c(positions[6,2]:positions[6,3]),]
Core.100x <- (sum(Corecoverage$V3 >=100)/nrow(Corecoverage))*100
Core.2x <- (sum(Corecoverage$V3 >1)/nrow(Corecoverage))*100
Core.mean <- round(mean(Corecoverage$V3),1)
Core.min  <- min(Corecoverage$V3)

# "ExLowRegion" Exclude End of RT to Start of PreC
ExLowRegcoverage <- coverage[c(1:positions[2,3],positions[5,2]:positions[8,3]),]
ExLowReg.100x <- (sum(ExLowRegcoverage$V3 >=100)/nrow(ExLowRegcoverage))*100
ExLowReg.2x <- (sum(ExLowRegcoverage$V3 >1)/nrow(ExLowRegcoverage))*100
ExLowReg.mean <- round(mean(ExLowRegcoverage$V3),1)
ExLowReg.min  <- min(ExLowRegcoverage$V3)

# "LowRegion" Only End of RT to Start of PreC
LowRegcoverage <- coverage[c(positions[2,3]:positions[5,2]),]
LowReg.100x <- (sum(LowRegcoverage$V3 >=100)/nrow(LowRegcoverage))*100
LowReg.2x <- (sum(LowRegcoverage$V3 >1)/nrow(LowRegcoverage))*100
LowReg.mean <- round(mean(LowRegcoverage$V3),1)
LowReg.min  <- min(LowRegcoverage$V3)




paste(WGS.100x, WGS.2x, WGS.mean, WGS.min, Pol.100x, Pol.2x, Pol.mean, Pol.min, Surface.100x, Surface.2x, Surface.mean, Surface.min, X.100x, X.2x, X.mean, X.min, PreC.100x, PreC.2x, PreC.mean, PreC.min, Core.100x, Core.2x, Core.mean, Core.min, ExLowReg.100x, ExLowReg.2x, ExLowReg.mean, ExLowReg.min, LowReg.100x, LowReg.2x, LowReg.mean, LowReg.min, collapse="\t")

CovReport <- data.frame(WGS.100x, WGS.2x, WGS.mean, WGS.min, Pol.100x, Pol.2x, Pol.mean, Pol.min, Surface.100x, Surface.2x, Surface.mean, Surface.min, X.100x, X.2x, X.mean, X.min, PreC.100x, PreC.2x, PreC.mean, PreC.min, Core.100x, Core.2x, Core.mean, Core.min, ExLowReg.100x, ExLowReg.2x, ExLowReg.mean, ExLowReg.min, LowReg.100x, LowReg.2x, LowReg.mean, LowReg.min)
paste(seqname,"F
