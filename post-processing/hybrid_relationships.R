library(stringr)
library(data.table)
library(plyr)
intable <- as.matrix(read.table("ubertree.tre",fill=TRUE))
suffixes <- as.numeric(as.matrix(read.table("suffix.txt")))

lenintable <- dim(intable)[1]
temptable <- NULL

for (j in 1:lenintable) {
temp <- unlist(strsplit(intable[j,2],"\\("))
for (k in 1:(length(temp))) {
temp1 <-  unlist(strsplit(temp[k],"\\)"))[1]
coloncount <- nchar(temp1) - nchar(gsub(":","",temp1,fixed=TRUE))
if (coloncount==2) {
temp2 <- unlist(strsplit(temp1,":"))
temp3 <- unlist(strsplit(temp2,","))
temp4 <- c(intable[j,1],temp3[1],temp3[3])
temptable <- rbind(temptable,temp4)
}
}
}

colnames(temptable) <- NULL
rownames(temptable) <- NULL

sumtable <- as.data.table(temptable)[,c(2,3) := list(pmin(V2, V3), pmax(V2, V3))]
sumtable$V2 <- substr(sumtable$V2,1, nchar(sumtable$V2)-suffixes)
sumtable$V3 <- substr(sumtable$V3,1, nchar(sumtable$V3)-suffixes)

out <- ddply(sumtable,.(V2,V3),nrow)

write.table(sumtable, "allele_combinations_by_locus.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(out, "allele_combinations.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()

