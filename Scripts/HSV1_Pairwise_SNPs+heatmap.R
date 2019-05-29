if(!require(Cairo)){
  install.packages("Cairo",repos="http://www.stats.bris.ac.uk/R/")
  library("Cairo")
}
if(!require(ggplot2)){
  install.packages("ggplot2",repos="http://www.stats.bris.ac.uk/R/")
  library(ggplot2)
}
###############
library(NMF)
library(RColorBrewer)
#########################

#basedropbox <- "C:/Users/Mat/Dropbox"
basedropbox <- "/Users/mb29/Dropbox"
#########


# Initially determined pairwise snps using MEGA - but could probably code this in for larger datasets
snps.columns <- read.csv(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Corededup+Szpara2014_v3s__SNPs_columns2.csv",sep=""), stringsAsFactors = F)

#snps.matrix <- read.matrix("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Corededup+Szpara2014_v3s__SNPs_Matrix.csv")
snps.columns$Dist <- as.numeric(snps.columns$Dist)

#snps.twice <- rbind(snps.columns, snps.columns[,c(2,1,3)])


p1 <- ggplot(snps.columns, aes(Species.1, Species.2))
p1 <- p1 +  geom_tile(aes(fill=Dist), color="white") +
  scale_fill_gradient(low="yellow",high="red", name="Pairwise\nSNPs") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(label = Dist), color = "black", size = 3) +
  coord_fixed()
p1  

Cairo(file="C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV+Szpara2014_Pairwise_SNPs_ggplot-heatmaptile.png" , width = 1200, height = 1100,type="png",dpi=900, units = "pt")
#aheatmap(as.matrix(d), color="YlOrRd", fontsize=14, cexRow=1, cexCol=1)
#aheatmap(d, color="YlOrRd", fontsize=14, cexRow=1, cexCol=1)
p1
dev.off()




### Now do for this study only (including duplicates)
#snps.columns.uk <- read.csv("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Core__SNPs_columns2.csv")
snps.columns.uk <- read.csv(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Core__SNPs_columns2.csv",sep=""), stringsAsFactors = F)


p2 <- ggplot(snps.columns.uk, aes(Species.1, Species.2))
p2 <- p2 +  geom_tile(aes(fill=Dist), color="white") +
  scale_fill_gradient(low="yellow",high="red", name="Pairwise\nSNPs") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(label = Dist), color = "black", size = 3) +
  coord_fixed()
p2  

Cairo(file="C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1-GoodCov_Core_Pairwise_SNPs_ggplot-heatmaptile.png" , width = 1200, height = 1100,type="png",dpi=900, units = "pt")
#aheatmap(as.matrix(d), color="YlOrRd", fontsize=14, cexRow=1, cexCol=1)
#aheatmap(d, color="YlOrRd", fontsize=14, cexRow=1, cexCol=1)
p2
dev.off()

###


### Now do for this study only (including duplicates)
#sites.columns.uk <- read.csv("C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Core__Sites_columns1.csv")
sites.columns.uk <- read.csv(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Core__Sites_columns1.csv",sep=""), stringsAsFactors = F)


p3 <- ggplot(sites.columns.uk, aes(Species.1, Species.2))
p3 <- p3 +  geom_tile(aes(fill=Dist), color="white") +
  scale_fill_gradient(low="yellow",high="red", name="Pairwise\nSNPs") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(label = Dist), color = "black", size = 1.5) +
  coord_fixed()
p3  

Cairo(file="C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1-GoodCov_Core_Pairwise_Sites_ggplot-heatmaptile.png" , width = 400, height = 300,type="png",dpi=600, units = "pt")
#aheatmap(as.matrix(d), color="YlOrRd", fontsize=14, cexRow=1, cexCol=1)
#aheatmap(d, color="YlOrRd", fontsize=14, cexRow=1, cexCol=1)
p3
dev.off()









#####
### Re-order according to patient (not randomly according to sample)
# Make double side
sites.columns.uk3 <- sites.columns.uk[,c(2,1,3)]
colnames(sites.columns.uk3) <- c("Species.1","Species.2","Dist")
sites.columns.uk2 <- rbind(sites.columns.uk, sites.columns.uk3)

### Re-order according to patient (not randomly according to sample)
patient.sample <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Sample+patient.tdl",sep=""), stringsAsFactors = F, header=T)
patient.sample$order <- c(1:nrow(patient.sample))
colnames(patient.sample) <- c("patient","seqname","order")
patient.sample$order <- factor(patient.sample$order,as.character(patient.sample$order))

#patient.sample$color <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44",  "#774411", "#AA7744", "#DDAA77", "#771122", "#DD7788")


sites.columns.uk2 <- merge(sites.columns.uk2, patient.sample, by.x="Species.1", by.y="seqname")
colnames(sites.columns.uk2)[5] <- "order.S1"

sites.columns.uk2 <- merge(sites.columns.uk2, patient.sample, by.x="Species.2", by.y="seqname") 
#colnames(sites.columns.uk2)[7] <- colnames(sites.columns.uk2)[5]
#sites.columns.uk2[,7]<- rev(sites.columns.uk2[,7])

#colnames(sites.columns.uk2)[11] <- "order.S2"
colnames(sites.columns.uk2)[7] <- "order.S2"

sites.columns.uk2 <- sites.columns.uk2[with(sites.columns.uk2,order(order.S1,order.S2)),]
#sites.columns.uk2$Species.1 <- factor(sites.columns.uk2$Species.1,as.character(sites.columns.uk2$Species.1))
#sites.columns.uk2$Species.2 <- factor(sites.columns.uk2$Species.2,as.character(sites.columns.uk2$Species.2))


sites.columns.uk2$Species.1 <- factor(sites.columns.uk2$Species.1, levels=(sites.columns.uk2$Species.1)[order(sites.columns.uk2$order.S1)])
sites.columns.uk2$Species.2 <- factor(sites.columns.uk2$Species.2, levels=(sites.columns.uk2$Species.2)[order(sites.columns.uk2$order.S2)])





p4 <- ggplot(sites.columns.uk2, aes(Species.1, Species.2))
p4 <- p4 +  geom_tile(aes(fill=Dist), color="white") +
  scale_fill_gradient(low="yellow",high="red", name="Pairwise\nSNPs") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(label = Dist), color = "black", size = 1.5) +
  coord_fixed()#geom_text(aes(label = SNPs), color = "black", size = 3) +
  theme(text = element_text(size=6)) +
  coord_fixed()
p4  

Cairo(file=paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1-GoodCov_Core_Pairwise_Sites_ggplot-heatmaptile-2ordered.png",sep="") , width = 400, height = 300,type="png",dpi=600, units = "pt")
p4
dev.off()




# Read in sample coloring file for axis labels
samples.colors <- read.table(paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1_Sample_colors2.tdl",sep=""), stringsAsFactors = F, header=F)
colnames(samples.colors) <- c("sample","color")
samples.colors <- data.frame(cbind(samples.colors,samples.colors$color[c(2:nrow(samples.colors),1)]),stringsAsFactors = F)
colnames(samples.colors) <- c("sample","colorx","colory")


#p4 + theme(axis.text.x = element_text(color=samples.colors$color))#, axis.text.y = element_text(samples.colors$color[c(2:nrow(samples.colors),1)]))
#p4 + theme( axis.text.y = element_text(samples.colors$color[c(2:nrow(samples.colors),1)]))

p4 <- p4 + theme(axis.text.x = element_text(color=samples.colors$colorx), axis.text.y = element_text(color=samples.colors$colorx))

Cairo(file=paste(basedropbox,"/UCL_Pathseek/HSV1/Ref_based/0.5_Snippy-cuttoff/SNP_stats/HSV1-GoodCov_Core_Pairwise_Sites_ggplot-heatmaptile-2ordered-coloured.png",sep="") , width = 400, height = 300,type="png",dpi=600, units = "pt")
p4
dev.off()

