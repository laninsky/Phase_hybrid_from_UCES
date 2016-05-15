library(stringr)
library(data.table)
library(plyr)
intable <- as.matrix(read.table("ubertree.tre",fill=TRUE))
suffixes <- as.numeric(as.matrix(read.table("suffix.txt")))

lenintable <- dim(intable)[1]
temptable <- NULL
temptablebysample <- NULL

for (j in 1:lenintable) {
temptemptable <- NULL
temp <- unlist(strsplit(intable[j,2],"\\("))
for (k in 1:(length(temp))) {
temp1 <-  unlist(strsplit(temp[k],"\\)"))[1]
coloncount <- nchar(temp1) - nchar(gsub(":","",temp1,fixed=TRUE))
if (coloncount==2) {
temp2 <- unlist(strsplit(temp1,":"))
temp3 <- unlist(strsplit(temp2,","))
temp4 <- c(intable[j,1],temp3[1],temp3[3])
temptemptable <- rbind(temptemptable,temp4)
temptable <- rbind(temptable,temp4)
}
}
colnames(temptemptable) <- NULL
rownames(temptemptable) <- NULL
sumtemptable <- as.data.table(temptemptable)[,c(2,3) := list(pmin(V2, V3), pmax(V2, V3))]
sumtemptable$V2 <- substr(sumtemptable$V2,1, nchar(sumtemptable$V2)-suffixes)
sumtemptable$V3 <- substr(sumtemptable$V3,1, nchar(sumtemptable$V3)-suffixes)
samples <- matrix(NA,ncol=1,nrow=((dim(sumtemptable)[1])*2))
samples[1:(dim(sumtemptable)[1]),1] <- sumtemptable$V2
samples[((dim(sumtemptable)[1])+1):(dim(samples)[1]),1] <- sumtemptable$V3
samplenames <- unique(samples)

temp1table <- matrix("NSSS",ncol=4,nrow=(dim(samplenames)[1]))
temp1table[,1] <- temptemptable[1,1]
temp1table[,2] <- samplenames

for (l in 1:(dim(temp1table)[1])) {
for (m in 1:(dim(temptemptable)[1])) {
if(grepl(temp1table[l,2],temptemptable[m,2],fixed=TRUE)) {
temp1table[l,3] <- temptemptable[m,3]
}
if(grepl(temp1table[l,2],temptemptable[m,3],fixed=TRUE)) {
temp1table[l,4] <- temptemptable[m,2]
}
}
}
temptablebysample <- rbind(temptablebysample,temp1table)
}

colnames(temptable) <- NULL
rownames(temptable) <- NULL
colnames(temptablebysample) <- NULL
rownames(temptablebysample) <- NULL

sumtable <- as.data.table(temptable)[,c(2,3) := list(pmin(V2, V3), pmax(V2, V3))]
sumtable$V2 <- substr(sumtable$V2,1, nchar(sumtable$V2)-suffixes)
sumtable$V3 <- substr(sumtable$V3,1, nchar(sumtable$V3)-suffixes)
out <- ddply(sumtable,.(V2,V3),nrow)
write.table(out, "allele_combinations.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

sumtable <- as.data.table(temptablebysample)[,c(3,4) := list(pmin(V3, V4), pmax(V3, V4))]
sumtable$V3 <- substr(sumtable$V3,1, nchar(sumtable$V3)-suffixes)
sumtable$V4 <- substr(sumtable$V4,1, nchar(sumtable$V4)-suffixes)
out <- ddply(sumtable,.(V3,V4),nrow)
write.table(sumtable, "allele_combinations_per_locus.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(out, "allele_combinations_per_locus_summary.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
