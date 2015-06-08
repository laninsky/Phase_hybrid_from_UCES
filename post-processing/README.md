#Post-processing#
Now we have our phased alleles for the hybrid tacked on to the fasta with the rest of our samples. In ONE_uce-xxx.fasta we have the first allele and the rest of our samples. In TWO_uce-xxx.fasta we have the second allele and the rest of our samples. You can go ahead and delete the rest of the stuff.

1) The first step is to re-align our hybrid alleles and the rest of the sequences, as indels were stripped out of the hybrid file so that gatk could eat the sequence for phasing. I am going to use MAFFT for this:
```
for i in *.fasta; do mv $i temp; mafft temp > $i; done;
```

2) Next, we need to convert our fasta alignments to phylip so we can run raxml. We'll use the phyluce wrappers for this (code quite liberarlly borrowed from Carl Oliveros - https://github.com/carloliveros/uce-scripts/blob/master/Species%20tree.md - thanks Carl!). Set the --alignments argument in the convert script to be where your fasta files are. You'll also want to select an appropriate outgroup for the 'phyluce_genetrees_run_raxml_genetrees.py' wrapper:
```
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignments uces --output phylip/ --input-format fasta --output-format phylip
python /public/uce/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_genetrees.py --input phylip --output phase_genetrees --outgroup k_med_d --cores 6 --quiet 2>&1 | tee kaloula_log/raxml_genetrees.txt
```

3) We then need to navigate to the folder with all our genetrees in it (in my case 'phase_genetrees'). We want to create a file similar to 'all-best-trees.tre', except with UCE names tied to each tree:
```
touch ubertree.tre
for i in `ls --color=never`; do if [ -d $i ]; then printname=`cat $i/RAxML_bestTree.best`; echo $i $printname >> ubertree.tre; fi; done 
```

4) Now we are going to use some R code to pull out the clades with just our hybrid in it. This assumes you've already installed the library package stringr, data.table and plyr e.g. install.packages("stringr"). It also assumes that all your names are in the format k_abc_x. You'll need to adjust the regular expressions if not. Make sure to change out "k_pix_e" for your actual hybrid's name...
```
library(stringr)
library(data.table)
library(plyr)
intable <- as.matrix(read.table("ubertree.tre"))
lenintable <- dim(intable)[1]
temptable <- matrix(NA, nrow=lenintable,ncol=3)

for (j in 1:lenintable) {
temptable[j,1] <- unlist(strsplit(intable[j,1],"_"))[2]
temptable[j,2] <- unlist(strsplit(intable[j,1],"_"))[1]
temptable[j,3] <- "NSSS"
temp <- unlist(strsplit((gsub(":[0-9]+\\.[0-9]+","",intable[j,2],fixed=FALSE)),"\\("))
lentemp <- length(temp)
for (i in 1:lentemp) {
if ((length(grep("k_pix_e",temp[i])))>0) {
if ((length(grep("k_[a-z]{3}_[a-z],k_pix_e\\))",temp[i])))>0) {
temptable[j,3] <- unlist(strsplit(temp[i],","))[1]
}
if ((length(grep("k_pix_e,k_[a-z]{3}_[a-z]\\))",temp[i])))>0) {
temptable[j,3] <- unlist(strsplit((gsub("k_pix_e,","",temp[i],fixed=TRUE)),"\\)"))[1]
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
```

5) All going well, there should now be a file in your working directory called "allele_combinations.txt". You can examine it and see what combinations of closest relatives for each allele were found across the gene-trees. "NSSS" is an abbreviation for "no single sister species" e.g. in this case the hybrid allele was sister to a clade containing multiple taxa. However, if your dataset is like mine, the fun doesn't stop there... potentially there other hybrids "messing" up your analysis that you want to remove before doing the analysis above. This code should prune your problem child and then run the same code above (NB: problem child here is k_cfc_r - replace this with your actual problem child, and don't forget to also replace the hybrid code - k_pix_e - with your sample of interest).
```
library(stringr)
library(data.table)
library(plyr)
intable <- as.matrix(read.table("ubertree.tre"))
lenintable <- dim(intable)[1]
temptable <- matrix(NA, nrow=lenintable,ncol=3)

for (j in 1:lenintable) {
temptable[j,1] <- unlist(strsplit(intable[j,1],"_"))[2]
temptable[j,2] <- unlist(strsplit(intable[j,1],"_"))[1]
temptable[j,3] <- "NSSS"
temp <- gsub(":[0-9]+\\.[0-9]+","",intable[j,2])
temp <- gsub("\\(k_pix_e,k_cfc_r\\)","k_pix_e",temp)
temp <- gsub("\\(k_cfc_r,k_pix_e\\)","k_pix_e",temp)


temp <- unlist(strsplit(temp,"\\("))
lentemp <- length(temp)
for (i in 1:lentemp) {
if ((length(grep("k_cfc_r",temp[i])))>0) {
if ((length(grep("k_[a-z]{3}_[a-z],k_cfc_r\\))",temp[i])))>0) {
temptable[j,3] <- unlist(strsplit(temp[i],","))[1]
}
if ((length(grep("k_cfc_r,k_[a-z]{3}_[a-z]\\))",temp[i])))>0) {
temptable[j,3] <- unlist(strsplit((gsub("k_cfc_r,","",temp[i],fixed=TRUE)),"\\)"))[1]
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
```



