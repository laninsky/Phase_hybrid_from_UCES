intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")
variables <- read.table("phasing_settings",header=FALSE,stringsAsFactors=FALSE,sep="\t")
namefile <- read.table("namefile",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]

no_taxa <- as.numeric(variables[3,1])*2
hybridname <- variables[4,1]

to_write <- matrix(NA,ncol=1,nrow=no_taxa)
to_write[1,1] <- intable[1,1]

to_write_title <- 2
sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- intable[j,1]
to_write_title <- to_write_title + 1
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[no_taxa,1] <- sequencepaste

to_write_temp <- matrix(NA,ncol=1,nrow=(no_taxa-2))
to_write_hybrid <-  matrix(NA,ncol=1,nrow=2)

rows <- dim(to_write)[1]
x <- 1
j <- 1

while (j < rows) {
if ((length(grep(hybridname,to_write[j,1])))>0) {
to_write_hybrid[1,1] <-  paste(">",namefile[1,1],sep="")
to_write_hybrid[2,1] <-  to_write[(j+1),1]
} else {
to_write_temp[x,1] <-  to_write[j,1]
to_write_temp[x+1,1] <-  to_write[(j+1),1]
x <- x + 2
}
j <- j + 2
}

write.table(to_write_temp, "temp.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(to_write_hybrid, "hybrid_ref.fa",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

q()
