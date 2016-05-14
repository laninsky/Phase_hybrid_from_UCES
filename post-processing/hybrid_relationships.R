library(stringr)
library(data.table)
library(plyr)
intable <- as.matrix(read.table("ubertree.tre",fill=TRUE))
suffixes <- as.numeric(as.matrix(read.table("suffix.txt")))

lenintable <- dim(intable)[1]
temptable <- NULL
temptemptable <- NULL

for (j in 1:lenintable) {
temp <- unlist(strsplit(intable[j,2],"\\("))
for (k in 1:(length(temp))) {
temp1 <-  unlist(strsplit(temp[k],"\\)"))[1]
coloncount <- nchar(temp1) - nchar(gsub(":","",temp1,fixed=TRUE))
if (coloncount==2) {
temp2 <- unlist(strsplit(temp1,":"))
temp3 <- unlist(strsplit(temp2,","))
temp4 <- c(intable[j,1],temp3[1],temp3[3])
temptemptable <- rbind(temptemptable,temp4)
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

temp1table <- matrix(NA,ncol=4,nrow=(dim(samplenames)[1]))
temp1table[,1] <- temptemptable[1,1]
temp1table[,2] <- samplenames





######## UP TO HERE RECODING ###########


rbind(,sumtemptable$V3)


}

colnames(temptable) <- NULL
rownames(temptable) <- NULL

sumtable <- as.data.table(temptable)[,c(2,3) := list(pmin(V2, V3), pmax(V2, V3))]

out <- ddply(sumtable,.(V2,V3),nrow)

write.table(sumtable, "allele_combinations_by_locus.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(out, "allele_combinations.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()

