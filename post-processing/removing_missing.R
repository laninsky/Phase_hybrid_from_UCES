filename <- as.matrix(read.table("name"))

temp <- readLines(filename)
notaxa <- as.numeric(unlist(strsplit(temp[1],"\\s+"))[2])
output <- matrix(0,ncol=1,nrow=notaxa)

for (i in 1:notaxa) {
dnasequence <- unlist(strsplit(temp[i+1],"\\s+"))
dnasequence <- paste(dnasequence[-1],temp[seq((2+i+notaxa),length(temp),(notaxa+1))],collapse="")
if ((grepl("A",dnasequence) | grepl("C",dnasequence) | grepl("G",dnasequence) | grepl("T",dnasequence))==TRUE) {
output[i,1] <- 1
}
}

firstrowoutput <- c(1,(which(output[,1]==1)+1))
finaloutput <- firstrowoutput
for (i in 1:((length(temp)/(notaxa+1))-1)) {
finaloutput <- c(finaloutput,(firstrowoutput+(i*(notaxa+1))))
}

finaloutput <- sort(finaloutput)
finaloutput <- temp[finaloutput]
firstrow <- unlist(strsplit(temp[1],"\\s+"))
firstrow[2] <- sum(output)
firstrow <- paste(firstrow,collapse=" ")
finaloutput[1] <- firstrow

outputname <- write.table(finaloutput,filename,quote=FALSE,row.names=FALSE,col.names=FALSE)
