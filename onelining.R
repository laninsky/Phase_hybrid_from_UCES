intable <- read.table("temp_alt.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]

to_write <- matrix(NA,ncol=1,nrow=2)
to_write[1,1] <- intable[1,1]

sequencepaste <- intable[2,1]

for (j in 3:rows) {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}

to_write[2,1] <- sequencepaste

write.table(to_write, "hybrid_ref2.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
