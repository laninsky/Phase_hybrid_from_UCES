library(stringr)
library(data.table)
library(plyr)
intable <- as.matrix(read.table("ubertree.tre",fill=TRUE))
species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE)

hybrid <- species[species[,2]=="hybrid",1]

lenintable <- dim(intable)[1]
temptable <- matrix(NA, nrow=lenintable,ncol=3)

hybfirst <-paste("*,",hybrid,"\\)",sep="")
hybsecond <-paste(hybrid,",.*","\\)",sep="")
hybreplace1 <- paste(",",hybrid,",", sep="")
hybreplace2 <- paste(hybrid,",", sep="")

for (j in 1:lenintable) {
temptable[j,1] <- unlist(strsplit(intable[j,1],"_"))[2]
temptable[j,2] <- unlist(strsplit(intable[j,1],"_"))[1]
temptable[j,3] <- "NSSS"
temp <- unlist(strsplit((gsub(":[0-9]+\\.[0-9]+","",intable[j,2],fixed=FALSE)),"\\("))
lentemp <- length(temp)
for (i in 1:lentemp) {
if ((length(grep(hybrid,temp[i])))>0) {
if ((length(grep(hybfirst,temp[i])))>0) {
temptable[j,3] <- unlist(strsplit(temp[i],","))[1]
}
if ((length(grep(hybsecond,temp[i])))>0) {
if ((length(grep(hybreplace1,temp[i])))>0) {
temptable[j,3] <- gsub("\\).*","",(unlist(strsplit((gsub(hybreplace2,",",temp[i],fixed=TRUE)),",,"))[2]))
} else {
temptable[j,3] <- unlist(strsplit((gsub(hybreplace2,"",temp[i],fixed=TRUE)),"\\)"))[1]
}
}
}
}
}

rm(intable)
rm(temp)
rm(lentemp)
rm(i)
rm(j)

temptable <- temptable[order(temptable[,1]),]
sumtable <- matrix(NA,ncol=3,nrow=(lenintable/2))

no_taxa <- dim(species)[1]

for (i in 1:lenintable) {
for (j in 1:no_taxa) {
if (temptable[i,3]==species[j,1]) {
if (!(species[j,2]=="hybrid")) {
temptable[i,3] <- species[j,2]
}
break
}
}
}

sumtable[1,1] <- temptable[1,1]
sumtable[1,2] <- temptable[1,3]
sumtable[1,3] <- temptable[2,3]

i <- 3
while (i < lenintable) {
j <- (i+1)/2
sumtable[j,1] <- temptable[i,1]
sumtable[j,2] <- temptable[i,3]
sumtable[j,3] <- temptable[i+1,3]
i <- i+2
}

sumtable2 <- as.data.table(sumtable)[,c("V2","V3") := list(pmin(V2, V3), pmax(V2, V3))]

out <- ddply(sumtable2,.(V2,V3),nrow)

write.table(out, "allele_combinations.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
q()

