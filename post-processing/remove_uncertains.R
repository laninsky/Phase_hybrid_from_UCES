intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")
species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]

no_final_taxa <- dim(species)[1] - sum(species[,2]=="uncertain")

to_write <- matrix(NA,ncol=1,nrow=(no_final_taxa*2))

uncertains <- species[species[,2]=="uncertain",1]
no_uncertains <- length(uncertains)

to_write_title <- 1

sequencepaste <- NULL

for (j in 1:rows) {
if ((length(grep(">",intable[j,1])))>0) {
if(!(any(uncertains %in% (gsub(">","",intable[j,1]))))) {
if(is.null(sequencepaste)) {
to_write[to_write_title,1] <- intable[j,1]
to_write_title <- to_write_title + 1
} else {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- intable[j,1]
to_write_title <- to_write_title + 1
}
}
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[(no_final_taxa*2),1] <- sequencepaste

write.table(to_write, "temp.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)
