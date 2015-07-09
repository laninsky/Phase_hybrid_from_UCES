intable1 <- read.table("hybrid_ref.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")
intable2 <- read.table("hybrid_ref2.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")
variables <- read.table("phasing_settings",header=FALSE,stringsAsFactors=FALSE,sep="\t")

hybridname <- paste(">",variables[4,1],sep="")

lengthtable <- dim(intable1)[1]

i <- 1

while (i < lengthtable) {
filename <- gsub(">","",intable1[i,1])
onefilename <- paste("ONE_",filename,sep="")
twofilename <- paste("TWO_",filename,sep="")

write.table(hybridname, onefilename, quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
write.table(intable1[i+1,1], onefilename, quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

write.table(hybridname, twofilename, quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
write.table(intable2[i+1,1], twofilename, quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

i <- i + 2
}

q()
