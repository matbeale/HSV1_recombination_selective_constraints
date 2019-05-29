

# dependencies
if(!require(Cairo)){
  install.packages("Cairo",repos="http://www.stats.bris.ac.uk/R/")
  library("Cairo")
}
if(!require(ggplot2)){
  install.packages("ggplot2",repos="http://www.stats.bris.ac.uk/R/")
  library(ggplot2)
}


# pull in gene positions file
gene.pos.file <- "C:/Users/Mat/Dropbox/UCL_Pathseek/HSV1/Ref_based/snippy/Variscan_selection/HSV1genes.bdf.modified"
gene.pos <- read.table(gene.pos.file, header =F, stringsAsFactors = F)
colnames(gene.pos) <- c("Start","End","Gene","Strand")

# generate a positional matrix for genes
#gene.matrix <- data.frame(c(1:152200), stringsAsFactors = F)
gene.matrix <- data.frame(seq(1,152200,50), stringsAsFactors = F)
colnames(gene.matrix) <- "Position"


# Create positional matrix with one column per gene entry - can take a long time.
testgene <- c(1:nrow(gene.pos))
for (line in testgene) {
  gene.matrix[line+1] <- sapply(1:nrow(gene.matrix), function(x) ifelse(gene.matrix[x,1] >=gene.pos[line,1] & gene.matrix[x,1] <= gene.pos[line,2] & gene.pos[line,4] == "+",  "+", ""))
  gene.matrix[line+1] <- sapply(1:nrow(gene.matrix), function(x) ifelse(gene.matrix[x,1] >=gene.pos[line,1] & gene.matrix[x,1] <= gene.pos[line,2] & gene.pos[line,4] == "-",  "-", gene.matrix[x,line+1]))
  colnames(gene.matrix)[line+1] <- gene.pos[line,3]
}

write.table(gene.matrix,file="HSV1_gene-matrix.txt", sep="\t", quote=T, col.names=T, row.names=F)
read.table("HSV1_gene-matrix.txt",header =T, stringsAsFactors = F)



gene.matrix.melted <- melt(gene.matrix, id.vars="Position")
gene.matrix.melted$value <- as.factor(gene.matrix.melted$value)


p2 <- ggplot(gene.matrix.melted, aes(x=Position, y=value, group=value, colour=value))
p2 <- p2 + geom_point(size=3, shape=15) + labs(x="Position", y="Strand") +
  scale_x_continuous(labels = scales::comma,breaks=seq(0, max(combined.variscan$Midpoint), round(max(combined.variscan$Midpoint) / 15,1 -nchar(ceiling(max(combined.variscan$Midpoint) / 15))))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(),strip.background=element_blank()) +
  scale_y_discrete() + ylim("-","+") + theme(legend.position = "none")
    #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
p2 



