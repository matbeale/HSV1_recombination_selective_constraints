
coveragepath <- "/Users/mb29/hsv1/coverage/"

list.file <- paste0(coveragepath,"list.txt")

coverage.list <- read.table(list.file,stringsAsFactors = F)$V1
totals <- 152222

coveragefile <- paste0(coveragepath,coverage.list[1])
coverage = read.table(coveragefile,header=F)  #Read in file
names(coverage) <- c("sample", "position","coverage")
#meancov.over1 <- sum(coverage$coverage)/nrow(coverage)


# "depth" file does not report 0 coverage positions - reintroduce these here
allpositions <- data.frame(c(1: totals[[1]]),stringsAsFactors = F)
#allpositions <- data.frame(c(1: 169175),stringsAsFactors = F)
names(allpositions) <- "position"
coverage <- merge(coverage,allpositions, by="position", all.y=T)
coverage$coverage <- as.numeric(coverage$coverage)
coverage$coverage <- as.numeric(ifelse(is.na(coverage$coverage),"0",coverage$coverage))
coverage$sample <- seqname


allpositions <- data.frame(c(1: totals[[1]]),stringsAsFactors = F)
names(allpositions) <- "position"
coverage.merge <- allpositions 
#coverage.merge <- NULL
for (seq in coverage.list){
  coveragefile <- paste0(coveragepath,seq)
  coverage <- read.table(coveragefile,header=F)
  names(coverage) <- c("sample", "position","coverage")
  coverage <- merge(coverage,allpositions, by="position", all.y=T)
  coverage <- data.frame(coverage$coverage,stringsAsFactors = F)
  colnames(coverage) <- gsub("\\.depth","",seq)
  coverage.merge <- cbind(coverage.merge,coverage)
}
coverage.merge$mean <- rowMeans(coverage.merge[,2:ncol(coverage.merge)], na.rm = FALSE, dims = 1)


#Poor Coverage regions
nocov <- coverage.merge[is.na(coverage.merge$mean),]

# 1-9183
# 71656-71788
# 117207-132551
# 143756-143836
# 145707-152222

#
