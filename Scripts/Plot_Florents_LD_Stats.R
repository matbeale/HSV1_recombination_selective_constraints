
library(ggplot2)
library(grid)
library(gridExtra)
library(Cairo)
library(reshape2)
library(gggenes)
library(ggforce)


#LD_r2 <- read.table("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/22-11-16/LDscan-all/LD_r2.LocalLD.minalfrq1.biallelicsites.max1gaps.tab")

LD_combined <- NULL
filelist <- c("LDscan-all","LDscan-CSF","LDscan-nonCSF")

#########
#basedropbox <- "C:/Users/Mat/Dropbox"
basedropbox <- "/Users/mb29/Dropbox"
#########

dataset <- "LDscan3000-all"

LD_r2_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/LD_r2.LocalLD-subsampled.minalfrq1.biallelicsites.max1gaps.tab",sep=""))
LD_fisher_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/LD_Fisher.LocalLD-subsampled.minalfrq1.biallelicsites.max1gaps.tab",sep=""))
#LD_physical_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/LDscan-all/physical-window2000_biallelicsites_density.tab",sep=""), header = T)
LD_physical_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/physical-window3000_biallelicsites_density.tab",sep=""), header = T)


LD_merge_all <- merge(LD_r2_all, LD_fisher_all, by="foci")
LD_merge_all <- merge(LD_merge_all, LD_physical_all, by="foci")

#LD_merge_all$sig <- ifelse(LD_merge_all$logcompldfisub >= 5,"sig","nonsig")
LD_merge_all$sig <- ifelse((LD_merge_all$logcompldfisub >= 5 & as.numeric(LD_merge_all$reportsnpdens) >=20) ,"sig","nonsig")

# edit 12/10/2017 - add in retention of significant SNP sites with <20 sites
LD_merge_all$sig <- sapply(1:nrow(LD_merge_all), function(x) ifelse((LD_merge_all$logcompldfisub[x] >= 5 & as.numeric(LD_merge_all$reportsnpdens[x]) >=10 & as.numeric(LD_merge_all$reportsnpdens[x]) <20) ,"weak",LD_merge_all$sig[x]))



##########################

# Bring in positional matrix data
gene.matrix <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1_gene-matrix.txt",sep=""),header =T, stringsAsFactors = F)
gene.matrix.melted <- melt(gene.matrix, id.vars="Position")
gene.matrix.melted$value <- as.factor(gene.matrix.melted$value)

# bring in gene names for text and determine midpoint in gene position for plotting
gene.pos.file <- paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified",sep="")
gene.pos <- read.table(gene.pos.file, header =F, stringsAsFactors = F, quote="")
colnames(gene.pos) <- c("Start","End","Gene","Strand")
gene.pos$Start <- as.numeric(gene.pos$Start)
gene.pos$End <- as.numeric(gene.pos$End)
gene.pos$midpoint <- as.numeric(gene.pos$Start) + ((as.numeric(gene.pos$End) - as.numeric(gene.pos$Start) )/2)

#keygenes <- c("UL9","UL16","UL19","UL42","UL43","UL51","US3","US4","US7","US8","US8A","US9","US10","US11")
keygenes <- c("UL8","UL17","UL18","UL41","UL42","UL43","UL46","UL50","UL51","UL52","US3","US7","US8","US8A","US9","US10","US10","US11","US12")

genetext <- gene.pos[gene.pos$Gene %in% keygenes,c("midpoint","Strand","Gene")]
genetext$text <- "genelabel"
genetext$text2 <- ifelse(genetext$Strand=="-", "Gene.rev","Gene.fwd")
#genetext <- genetext[,c(1,4,3)]
genetext <- genetext[,c(1,5,3)]
colnames(genetext) <- c("Position","variable","value")
genetext$Position <- as.numeric(genetext$Position)
gene.matrix.melted <- rbind(gene.matrix.melted,genetext)
gene.matrix.melted$Position <- as.numeric(gene.matrix.melted$Position)

#  "US8","US8A","US9","US10","US11"s
gene.matrix.melted$jit <- with(gene.matrix.melted, ifelse((value  =="US8A" | value == "US10"), 1, 0.01))


# plot gene position layout 'p3'
p3 <- ggplot(gene.matrix.melted)
p3 <- p3 + geom_point(size=5, data=gene.matrix.melted[!(grepl("Gene",gene.matrix.melted$variable)),], aes(x=Position,y=value, group=value, colour=value,shape=value)) + labs(x="Genome Position (bp)", y="Strand") +
  #scale_x_continuous(limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  scale_x_continuous(expand=c(0,0),limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  theme_bw() + theme(legend.position = "none") +
  #scale_y_discrete() + ylim("-","+","Gene") +
  #scale_y_discrete() + ylim("-","+","genepos","geneneg") +
  scale_y_discrete() + ylim("Gene.rev","-","+","Gene.fwd") +
  #scale_shape_manual(values=c(0,"\u25BA","\u25C4")) # Try triangles - doesn't seem to print with Cairo as non-standard unicode
  scale_shape_manual(values=c(0,60, 62,0)) # Use '<' and '>' for gene positions
# add in text to plot
p3 <- p3 + geom_text(data=gene.matrix.melted[(grepl("Gene",gene.matrix.melted$variable)),], aes(x=Position, y=variable, label=value), angle = -90, size=2, alpha=1/2,
                     position=position_jitter(width=(gene.matrix.melted[(grepl("Gene",gene.matrix.melted$variable)),]$jit)/1,height=(gene.matrix.melted[(grepl("Gene",gene.matrix.melted$variable)),]$jit)*0.5))
p3


##### # Alternate way of plotting gene pos (p3) using new gggenes package
#library("gggenes")
gene.pos$Start.dir <- ifelse(gene.pos$Strand=='+', gene.pos$Start, gene.pos$End)
gene.pos$End.dir <- ifelse(gene.pos$Strand=='+', gene.pos$End, gene.pos$Start)
#gene.pos$keygenes <- 

genetext = NULL
genetext <- data.frame(gene.pos[gene.pos$Gene %in% keygenes,c("midpoint","Strand","Gene","Start.dir","End.dir")],stringsAsFactors = F)
genetext$jit <- with(genetext, ifelse((Gene  =="US8A" | Gene == "US10"), 0.5, 0.01))

#genetext <- gene.pos[gene.pos$Gene %in% keygenes,c("midpoint","Strand","Gene")]
  
#p3 <- ggplot(gene.pos,aes(xmin = Start.dir, xmax = End.dir, y = "molecule", fill = Gene)) +
#p3 <- ggplot(gene.pos,aes(xmin = Start.dir, xmax = End.dir, y = "molecule")) +  
p3 <- ggplot(gene.pos,aes(xmin = Start.dir, xmax = End.dir, y = "molecule", fill = Strand)) +
    geom_gene_arrow() + 
  scale_x_continuous(expand=c(0,0),limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) + 
  theme_bw() + theme(legend.position="none") +
  theme( axis.text.y = element_blank(), axis.ticks.y= element_blank()) + labs(y="Genes")
p3 <- p3 + geom_text(data=genetext, aes(x=midpoint, y="Genes", label=Gene),angle = -90, size=2, alpha=1/2,position=position_jitter(height = genetext$jit))
p3
###########################


# Plot LD squared by position (single run, e.g. 'All')
pld <- ggplot(LD_merge_all, aes(foci,meanldrsub, colour=sig))
pld <- pld + geom_point(alpha=1/4, size=2) + theme_bw() +
  scale_x_continuous(expand=c(0,0),limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  xlab("Genome Position (3000bp windows)") + ylab(expression(paste("LD Strength (r" ^{"2"}, ")", sep=""))) +
  scale_colour_manual(values=c("grey45","red","Blue3")) + theme(legend.position = "none") 
#  scale_colour_manual(values=c("grey45","red")) + theme(legend.position = "none") 
#pld
p4 <- ggplot(gene.pos,aes(xmin = Start.dir, xmax = End.dir, y = "molecule", fill = Strand)) +
  geom_gene_arrow() + 
  scale_x_continuous(expand=c(0,0),limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) + 
  theme_bw() + theme(legend.position="none") +
  theme( axis.text.y = element_blank(), axis.ticks.y= element_blank()) + labs(y="Genes") +
  geom_text(data=gene.pos, aes(x=midpoint, y="Genes", label=Gene),angle = -90, size=2, alpha=1/2,position=position_jitter(height = genetext$jit))

gA <- ggplotGrob(pld + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x= element_blank()))
gB <- ggplotGrob(p4 +xlab("Genome Position (3000bp windows)")) #+geom_vline(xintercept=143317, colour="blue3", alpha=1/3) + geom_vline(xintercept=93112, colour="blue3", alpha=1/3))
gB$widths <- gA$widths 
grid.newpage()
grid.arrange(gA, gB, nrow = 2, heights=c(8,2)) # for gggenes plotting


#pldfish <- ggplot(LD_fisher, aes(foci,logcompldfisub))
#pldfish + geom_point(alpha=1/5) + theme_bw() +
#  scale_x_continuous(limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
#  xlab("Genome Position (2000bp windows)") + ylab("LD Significance (-log10(p))")






############ Loop all 3 together
LD_combined <- NULL
#filelist <- c("LDscan-all","LDscan-CSF","LDscan-nonCSF")
filelist <- c("LDscan3000-all","LDscan3000-CSF","LDscan3000-nonCSF")
for (dataset in filelist){

LD_r2_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/LD_r2.LocalLD-subsampled.minalfrq1.biallelicsites.max1gaps.tab", sep=""))
LD_fisher_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/LD_Fisher.LocalLD-subsampled.minalfrq1.biallelicsites.max1gaps.tab", sep=""))
#LD_physical_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/physical-window2000_biallelicsites_density.tab", sep=""), header = T)
LD_physical_all <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/",dataset,"/physical-window3000_biallelicsites_density.tab", sep=""), header = T)

LD_merge_all <- merge(LD_r2_all, LD_fisher_all, by="foci")
LD_merge_all <- merge(LD_merge_all, LD_physical_all, by="foci")
#LD_merge_all <- merge(LD_r2_all, LD_physical_all, by="foci")

LD_merge_all$sig <- ifelse(LD_merge_all$logcompldfisub >= 5,"sig","nonsig")
LD_merge_all$sig <- ifelse((LD_merge_all$logcompldfisub >= 5 & as.numeric(LD_merge_all$reportsnpdens) >=20) ,"sig","nonsig")
#LD_merge_all$sig <- sapply(1:nrow(LD_merge_all), function (x) ifelse(as.numeric(LD_merge_all$reportsnpdens[x]) <20 ,"nonsig",LD_merge_all$sig[x]))

###
# edit 12/10/2017 - add in retention of significant SNP sites with <20 sites
LD_merge_all$sig <- sapply(1:nrow(LD_merge_all), function(x) ifelse((LD_merge_all$logcompldfisub[x] >= 5 & as.numeric(LD_merge_all$reportsnpdens[x]) >=10 & as.numeric(LD_merge_all$reportsnpdens[x]) <20) ,"weak",LD_merge_all$sig[x]))
###

#######
# Deal with low snp coverage bins for plotting

#LD_merge_all$sig <- sapply(1:nrow(LD_merge_all), function (x) ifelse(as.numeric(LD_merge_all$reportsnpdens[x]) <20 ,"lowcov",LD_merge_all$sig[x]))
#LD_merge_all$cov <-  ifelse(as.numeric(LD_merge_all$reportsnpdens) <20  ,"low","good")
LD_merge_all$meanldrsub <-  sapply(1:nrow(LD_merge_all), function(x) ifelse(as.numeric(LD_merge_all$reportsnpdens[x]) <10  ,"",LD_merge_all$meanldrsub[x]))
LD_merge_all$meanldrsub <- as.numeric(LD_merge_all$meanldrsub)
LD_merge_all$sig <-  sapply(1:nrow(LD_merge_all), function(x) ifelse(as.numeric(LD_merge_all$reportsnpdens[x]) <10  ,"",LD_merge_all$sig[x]))

#LD_merge_all$sig <- ifelse(as.numeric(LD_merge_all$reportsnpdens) >=20 ,"sig","nonsig")
#LD_merge_all$sig <- ifelse((as.numeric(LD_merge_all$meanldrsub) >= 0.15 & as.numeric(LD_merge_all$reportsnpdens) >=20) ,"sig","nonsig")

#LD_merge_all$sig <- ifelse((as.numeric(LD_merge_all$meanldrsub) >= 0.15 & as.numeric(LD_merge_all$reportsnpdens) >=20) ,"sig",LD_merge_all$sig)
############

LD_merge_all$dataset <- gsub("LDscan-","",dataset)

LD_combined <- rbind(LD_combined,LD_merge_all)
}

LD_combined$dataset <- gsub("all", "All",LD_combined$dataset)

# Rename Groups for plot
#LD_combined$dataset <- gsub("^.+nonCSF", "Vesicle",LD_combined$dataset)
LD_combined$dataset <- gsub("^.+nonCSF", "Swab",LD_combined$dataset)
LD_combined$dataset <- gsub("^.+CSF", "CSF",LD_combined$dataset)
LD_combined$dataset <- gsub("^.+All", "All",LD_combined$dataset)


pld.loop <- ggplot(LD_combined, aes(foci,meanldrsub, colour=sig))
pld.loop <- pld.loop + geom_point(alpha=1/10, size=1) + theme_bw() +
  #scale_x_continuous(limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  scale_x_continuous(expand=c(0,0),limits = c(0, 152222),labels = scales::comma,breaks=seq(0, 152222, round(152222 / 15,1 -nchar(ceiling(152222 / 15))))) +
  #xlab("Genome Position (2000bp windows)") + ylab(expression(paste("LD Strength (r" ^{"2"}, ")", sep=""))) +
  xlab("Genome Position (3000bp windows)") + ylab(expression(paste("LD Strength (r" ^{"2"}, ")", sep=""))) +
  #ylim(0,0.8) +
  #ylim(0,0.7) +
  scale_y_continuous(expand=c(0,0),limits = c(0, 0.7)) +
  #scale_colour_manual(values=c("black","red")) + theme(legend.position = "none") +
  #scale_colour_manual(values=c("white","black","red")) + theme(legend.position = "none") +
  scale_colour_manual(values=c("white","black","red","blue3")) + theme(legend.position = "none") +
  #geom_tile(aes(x=foci, y=0.5, fill=cov)) +
  facet_grid(dataset ~ ., scales="free", space="free", drop=T, margins=F) +
  theme(strip.background = element_rect(fill="white"))
pld.loop

# add a vertical line to illustrate US9 start position (143317) or US8 end position (142899)
# start of UL42 =93112
#pld.loop <- pld.loop + geom_vline(xintercept=142899, colour="blue3", alpha=1/3)
#p3 <- p3 + geom_vline(xintercept=142899, colour="blue3", alpha=1/3)



# Use ggplotGrob to alter the widths so they line up, then plot with grid.arrange
gA <- ggplotGrob(pld.loop + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x= element_blank())+geom_vline(xintercept=142899, colour="blue3", alpha=1/3)+ geom_vline(xintercept=93112, colour="blue3", alpha=1/3))
#gB <- ggplotGrob(p3 +xlab("Genome Position (2000bp windows)")+geom_vline(xintercept=143317, colour="blue3", alpha=1/3) + geom_vline(xintercept=93112, colour="blue3", alpha=1/3))

gB <- ggplotGrob(p3 + geom_vline(xintercept=142899, colour="blue3", alpha=1/3)+ geom_vline(xintercept=93112, colour="blue3", alpha=1/3) +xlab("Genome Position (3000bp windows)")) #+geom_vline(xintercept=143317, colour="blue3", alpha=1/3) + geom_vline(xintercept=93112, colour="blue3", alpha=1/3))
#gB <- ggplotGrob(p3 +xlab("Genome Position (3000bp windows)"))
gB$widths <- gA$widths 
grid.newpage()



#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_r2-by-FisherSig.LocalLD-subsampled_+genes.png",sep=""), width = 800, height = 550,type="png",dpi=600, units = "pt")
#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_r2-by-FisherSig.LocalLD-subsampled_+genes.png",sep=""), width = 600, height = 420,type="png",dpi=600, pointsize = 11*600/72, units = "pt", family="Symbola")
#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_r2-by-FisherSig.LocalLD-subsampled_+genes.png",sep=""), width = 600, height = 420,type="png",dpi=600, pointsize = 10*600/72, units = "pt", family="Symbola")

#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_r2-by-FisherSig.LocalLD-subsampled_+genes.png",sep=""), width = 650, height = 420,type="png",dpi=600, pointsize = 10*600/72, units = "pt", family="Symbola")

Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_r2-by-FisherSig.LocalLD-subsampled_3000win+10to20-blue_+genes.png",sep=""), width = 650, height = 420,type="png",dpi=600, pointsize = 10*600/72, units = "pt", family="Symbola")

#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_3000bp_r2-by-FisherSig.LocalLD-subsampled_+genes.png",sep=""), width = 800, height = 550,type="png",dpi=600, units = "pt")
#grid.arrange(gA, gB, nrow = 2, heights=c(8,3))
#grid.arrange(gA, gB, nrow = 2, heights=c(6,2)) # for manually plotted

grid.arrange(gA, gB, nrow = 2, heights=c(8,2)) # for gggenes plotting

dev.off()







############
# Explore differences between samples
#
#sig_LD_combined <- subset(LD_combined, sig=="sig") # ? best to filter beforehand, or apply new filters after merge for direct comparison...?
sig_LD_combined <- LD_combined

sig_LD_combined.all <- subset(sig_LD_combined, dataset=="All")
sig_LD_combined.all$genemappings <- gsub("^, ","",sapply(1:nrow(sig_LD_combined.all), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(sig_LD_combined.all$foci[y] >= gene.pos$Start[x] && sig_LD_combined.all$foci[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)

sig_LD_combined.CSF <- subset(sig_LD_combined, dataset=="CSF")
sig_LD_combined.CSF$genemappings <- gsub("^, ","",sapply(1:nrow(sig_LD_combined.CSF), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(sig_LD_combined.CSF$foci[y] >= gene.pos$Start[x] && sig_LD_combined.CSF$foci[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)

#sig_LD_combined.nCSF <- subset(sig_LD_combined, dataset=="nonCSF")
#sig_LD_combined.nCSF <- subset(sig_LD_combined, dataset=="Vesicle")
sig_LD_combined.nCSF <- subset(sig_LD_combined, dataset=="Swab")
sig_LD_combined.nCSF$genemappings <- gsub("^, ","",sapply(1:nrow(sig_LD_combined.nCSF), function(y) toString(unique(sapply(1:nrow(gene.pos), function(x) ifelse(sig_LD_combined.nCSF$foci[y] >= gene.pos$Start[x] && sig_LD_combined.nCSF$foci[y] <= gene.pos$End[x], gene.pos$Gene[x],''))))),perl=T)


# Determine mean/median (baseline) LD for each dataset to allow normalisation
LD_all_median <- median(as.numeric(sig_LD_combined.all$meanldrsub), na.rm=T)
LD_CSF_median <- median(as.numeric(sig_LD_combined.CSF$meanldrsub), na.rm=T)
LD_nCSF_median <- median(as.numeric(sig_LD_combined.nCSF$meanldrsub), na.rm=T)

sig_LD_combined.CSF$normalisedldr.CSF <- sig_LD_combined.CSF$meanldrsub/LD_CSF_median
sig_LD_combined.nCSF$normalisedldr.nCSF <- sig_LD_combined.nCSF$meanldrsub/LD_nCSF_median

# Try and combine datasets (CSF/nonCSF)
sig_LD_CSF<- sig_LD_combined.CSF[,c("foci","meanldrsub","logcompldfisub","genemappings","reportsnpdens","normalisedldr.CSF")]
colnames(sig_LD_CSF) <- c("foci","CSF.meanldrsub","CSF.logcompldfisub", "genemappings","CSF.snpdens","normalisedldr.CSF")
sig_LD_nCSF<- sig_LD_combined.nCSF[,c("foci","meanldrsub","logcompldfisub","genemappings","reportsnpdens","normalisedldr.nCSF")]
colnames(sig_LD_nCSF) <- c("foci","nCSF.meanldrsub","nCSF.logcompldfisub", "genemappings","nCSF.snpdens","normalisedldr.nCSF")

sig_LD_merge <- merge(sig_LD_CSF,sig_LD_nCSF,by="foci", all.x=T, all.y=T)

sig_LD_merge <- sig_LD_merge[!(is.na(sig_LD_merge$CSF.meanldrsub) & is.na(sig_LD_merge$nCSF.meanldrsub)),] # remove sites absent from both datasets
sig_LD_merge <- sig_LD_merge[((sig_LD_merge$CSF.snpdens>=10) & (sig_LD_merge$nCSF.snpdens>=10)),] # only keep sites where at least 10 sites/window (min criteria for p value is 20, but want to include sites with less in case they are high anyway for qc puposes)
sig_LD_merge <- sig_LD_merge[((sig_LD_merge$CSF.logcompldfisub>=5) | (sig_LD_merge$nCSF.logcompldfisub>=5)),] # only keep sites where at least one is significant

# Determine the differences between the CSF and non-CSF
#sig_LD_merge$difference <- sapply(1:nrow(sig_LD_merge), function(x) diff(as.numeric(sig_LD_merge[x,c("CSF.meanldrsub","nCSF.meanldrsub")])))
sig_LD_merge$difference <- sapply(1:nrow(sig_LD_merge), function(x) diff(as.numeric(sig_LD_merge[x,c("normalisedldr.CSF","normalisedldr.nCSF")])))
# Relabel the 'unmapped'blank hits as "No Gene"
sig_LD_merge$genemappings <- sapply(1:nrow(sig_LD_merge), function(x) ifelse(sig_LD_merge$genemappings.y[x]=="","Non Coding",sig_LD_merge$genemappings.y[x]))

#############

# Force order of plot
#gene.order <- data.frame(rev(c("UL9","UL16","UL19","UL42","UL43","UL50","UL51","US3","US4","US7","US8","US8, US8A","US8A","US9","US10","US10, US11","US11","Non Coding")), stringsAsFactors= T)
gene.order <- data.frame(rev(c("UL8","UL17","UL18","UL41","UL42","UL43","UL46","UL50","UL51","UL52","US3","US7","US8","US8A","US9","US10","US10, US11","US11","US12","Non Coding")), stringsAsFactors= T)

gene.order$order <- factor(c(1:nrow(gene.order)))
colnames(gene.order) <- c("genelab","order")

gene.order$order2 <- factor(gene.order$order,levels = gene.order$order[order(gene.order$order)])
#gene.order$order <- factor(gene.order$order,as.character(gene.order$order))
#gene.order$order <- factor(as.character(gene.order$order))

sig_LD_merge2 <- merge(sig_LD_merge, gene.order, by.x="genemappings", by.y="genelab")
sig_LD_merge2 <- sig_LD_merge2[with(sig_LD_merge2,order(order)),]
sig_LD_merge2$genemappings <- as.character(sig_LD_merge2$genemappings)

sig_LD_merge2$genemappings <- factor(sig_LD_merge2$genemappings, levels = gene.order$genelab[order(gene.order$order)])

#########

###
# Plot comparison analysis figure

#library(ggforce)
pmappings <- ggplot(sig_LD_merge2, aes(genemappings, difference, colour=genemappings)) 
pmappings <- pmappings + geom_sina(alpha=1/3) + theme_bw() +
  #ylab("\n\nNormalised LD difference\n(Sig hits; nCSF-CSF)")+  #ylab(expression(paste("\nNormalised LD r" ^{"2"}, "difference\n(nCSF-CSF)", sep="")))
  ylab("\n\nNormalised LD difference\n(Sig hits; SWAB-CSF)")+  #ylab(expression(paste("\nNormalised LD r" ^{"2"}, "difference\n(nCSF-CSF)", sep=""))) 
  xlab("Gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5)) +
  theme(legend.position = "none") + geom_hline(yintercept=0) 
pmappings + coord_flip() 

#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_CSF-vs_nCSF_gene-comparison.png",sep=""), width = 800, height = 200,type="png",dpi=600, units = "pt")
#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_CSF-vs_nCSF_gene-comparison.png",sep=""), width = 200, height = 400,type="png",dpi=600,pointsize = 11*600/72, units = "pt")
Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_CSF-vs_Vesicle_gene-comparison__3000win_16-10-17.png",sep=""), width = 200, height = 400,type="png",dpi=600,pointsize = 11*600/72, units = "pt")
pmappings + coord_flip() + scale_x_discrete(position = "top")
dev.off()
#############


# Combine this below the main plot?? - no, doesn't look good - better to have it rotated on the side...
#gC <- ggplotGrob(pmappings)
#gC$widths <- gA$widths 
#Cairo(file= paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/Florent__LD/19-01-17/HSV1_combined_LD_r2-by-FisherSig.LocalLD-subsampled_+genes_+_Gene-differences.png",sep=""), width = 800, height = 800,type="png",dpi=600, units = "pt")
#grid.arrange(gA, gB, gC, nrow = 3, heights=c(5,1,2))
#dev.off()
