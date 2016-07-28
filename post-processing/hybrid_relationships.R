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
temp4 <- c(intable[j,1],temp3[1],temp3[3], NA)
temptemptable <- rbind(temptemptable,temp4)
}
}

for (i in 1:(dim(temptemptable)[1])) {
if(substr(temptemptable[i,2],1,(nchar(temptemptable[i,2])-suffixes))==substr(temptemptable[i,3],1,(nchar(temptemptable[i,3])-suffixes))) {
findpattern <- paste("\\(",temptemptable[i,2],":0.[0-9]+,",temptemptable[i,3],":0.[0-9]+\\)",sep="")
replacepattern <- paste((substr(temptemptable[i,2],1,(nchar(temptemptable[i,2])-suffixes))),paste(rep("x",suffixes),collapse=""),sep="")
intable[j,2] <- gsub(findpattern,replacepattern,intable[j,2])
}
}

temp <- unlist(strsplit(intable[j,2],"\\("))
for (k in 1:(length(temp))) {
temp1 <-  unlist(strsplit(temp[k],"\\)"))[1]
coloncount <- nchar(temp1) - nchar(gsub(":","",temp1,fixed=TRUE))
if (coloncount==2) {
temp2 <- unlist(strsplit(temp1,":"))
temp3 <- unlist(strsplit(temp2,","))
temp1lines <- which((substr(temp3[1],1,(nchar(temp3[1])-suffixes)))==substr(temptemptable[,2],1,(nchar(temptemptable[,2])-suffixes)))
temp3lines <- which((substr(temp3[3],1,(nchar(temp3[3])-suffixes)))==substr(temptemptable[,2],1,(nchar(temptemptable[,2])-suffixes)))
if(length(temp1lines)>0 && (substr(temptemptable[temp1lines,2],1,(nchar(temptemptable[temp1lines,2])-suffixes))==substr(temptemptable[temp1lines,3],1,(nchar(temptemptable[temp1lines,3])-suffixes)))) {
temptemptable[temp1lines,4] <- temp3[3]
}
if(length(temp3lines)>0 && (substr(temptemptable[temp3lines,2],1,(nchar(temptemptable[temp3lines,2])-suffixes))==substr(temptemptable[temp3lines,3],1,(nchar(temptemptable[temp3lines,3])-suffixes)))) {
temptemptable[temp3lines,4] <- temp3[1]
}
if(!(temp3[1] %in% temptemptable[,4])) {
if(!((temp3[1] %in% temptemptable[,2]) || (temp3[1] %in% temptemptable[,3]))) {
temp4 <- c(intable[j,1],temp3[1],temp3[3], NA)
temptemptable <- rbind(temptemptable,temp4)
}
}
if(!(temp3[3] %in% temptemptable[,4])) {
if(!((temp3[3] %in% temptemptable[,2]) || (temp3[3] %in% temptemptable[,3]))) {
temp4 <- c(intable[j,1],temp3[1],temp3[3], NA)
temptemptable <- rbind(temptemptable,temp4)
}
}
}
}

colnames(temptemptable) <- NULL
rownames(temptemptable) <- NULL
sumtemptable <- as.data.table(temptemptable)[,c(2,3) := list(pmin(V2, V3), pmax(V2, V3))]
sumtemptable$V2 <- substr(sumtemptable$V2,1, nchar(sumtemptable$V2)-suffixes)
sumtemptable$V3 <- substr(sumtemptable$V3,1, nchar(sumtemptable$V3)-suffixes)
sumtemptable$V4 <- substr(sumtemptable$V4,1, nchar(sumtemptable$V4)-suffixes)
samples <- matrix(NA,ncol=1,nrow=((dim(sumtemptable)[1])*2))
samples[1:(dim(sumtemptable)[1]),1] <- sumtemptable$V2
samples[((dim(sumtemptable)[1])+1):(dim(samples)[1]),1] <- sumtemptable$V3
samplenames <- unique(samples)

temp1table <- matrix("NSSS",ncol=5,nrow=(dim(samplenames)[1]))
temp1table[,1] <- temptemptable[1,1]
temp1table[,2] <- samplenames
temp1table[,5] <- NA

for (l in 1:(dim(temp1table)[1])) {
x <- 1
for (m in 1:(dim(temptemptable)[1])) {
if (temp1table[l,2]==substr(temptemptable[m,3],1,(nchar(temptemptable[m,3])-suffixes)) && substr(temptemptable[m,2],1,(nchar(temptemptable[m,2])-suffixes))==substr(temptemptable[m,3],1,(nchar(temptemptable[m,3])-suffixes))) {
temp1table[l,3:5] <- temptemptable[m,2:4]
break
} else {
if(temp1table[l,2]==(substr(temptemptable[m,2],1,(nchar(temptemptable[m,2])-suffixes)))) {
temp1table[l,(x+2)] <- temptemptable[m,3]
x <- x+1
}
if(temp1table[l,2]==(substr(temptemptable[m,3],1,(nchar(temptemptable[m,3])-suffixes)))) {
temp1table[l,(x+2)] <- temptemptable[m,2]
x <- x+1
}
}
}
}

temptablebysample <- rbind(temptablebysample,temp1table)
temptable <- rbind(temptable,temptemptable)
}

colnames(temptablebysample) <- NULL
rownames(temptablebysample) <- NULL
title <- c("locus","allele1","allele2","next_closest_allele")
temptable <- rbind(title,temptable)
write.table(temptable, "raw_allele_combinations.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

sumtable <- as.data.table(temptablebysample)[,c(3,4) := list(pmin(V3, V4), pmax(V3, V4))]
sumtable$V3 <- substr(sumtable$V3,1, nchar(sumtable$V3)-suffixes)
sumtable$V4 <- substr(sumtable$V4,1, nchar(sumtable$V4)-suffixes)
sumtable$V5 <- substr(sumtable$V5,1, nchar(sumtable$V5)-suffixes)
out <- ddply(sumtable,.(V2,V3,V4,V5),nrow)
title <- c("sample","sister_species1","sister_species2","next_closest_species","count")
out <- rbind(title,out)
write.table(out, "summmarized_allele_combinations.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

title <- c("locus","sample","allele1","allele2","next_closest_allele")
temptablebysample <- rbind(title,temptablebysample)
write.table(temptablebysample, "allele_combinations_by_locus.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
