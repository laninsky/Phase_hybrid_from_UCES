library(stringr)
intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")
rows <- dim(intable)[1]

locus_count <- 1
tempfile <- NULL

for (j in 1:rows) {
if ((length(grep("//",intable[j,1])))>0) {
filename <- paste(locus_count,".fasta",sep="")
write.table(tempfile, filename,quote=FALSE, col.names=FALSE,row.names=FALSE)
locus_count <- locus_count + 1
tempfile <- NULL
} else {
tempy <- NULL
tempy <- rbind(unlist(strsplit(intable[j,1],"[\\s ]+",fixed=FALSE))[1],unlist(strsplit(intable[j,1],"[\\s ]+",fixed=FALSE))[2])
tempfile <- rbind(tempfile,tempy)
}
}
q()
