#Post-processing#
Now we have our phased alleles for the hybrid tacked on to the fasta with the rest of our samples. In ONE_uce-xxx.fasta we have the first allele and the rest of our samples. In TWO_uce-xxx.fasta we have the second allele and the rest of our samples.

1) The first step is to copy all of these files to a new folder, just in case we stuff up at any of the downstream steps, so that we don't have to repeat everything
```
mkdir fasta
cp *.fasta fasta
```

2) Now we are going to re-align our hybrid alleles and the rest of the sequences, as indels were stripped out of the hybrid file so that gatk could eat the sequence for phasing. I am going to use MAFFT for this, but you could pick your favorite aligher:
```
cd fasta
for i in *.fasta; do mv $i temp; mafft temp > $i; done;
rm temp
```

3) Now, you might be in the situation where you have multiple hybrids/samples of uncertain taxonomic status in your dataset. These are going to muddle things up a little, because if the hybrid you have phased maps back to one of these 'uncertains', that doesn't give you a lot of information on the putative parental species of the hybrid you are interested in. So... at this point we can run remove_uncertains.R before we move on to the raxml stage. To do this, you'll need to set up the file 'species_assignments'. You might as well do this now, anyway, because you are going to need it later on. This file has the sample names in the first column (tab delimited), and the species/taxon designations you'd like to use in the second column. Make sure to put 'hybrid' for your putative hybrid, and 'uncertain' for any samples you think should be removed from your dataset e.g.
```
k_bal_r baleata
k_cfc_r uncertain
k_war_r k_war_r
k_con_r conjuncta
k_wal_r uncertain
k_pixme hybrid
k_mer_r meridionalis
k_koc_r kocakii
k_bal_r baleata
```

Upload species_assignments and remove_uncertains.R to the fasta directory which has your aligned fasta files in it (e.g. working_dir/fasta) and run it by:
```
for i in `ls *.fasta`;
do mv $i temp;
Rscript remove_uncertains.R
mv temp.fa $i;
rm -rf temp;
done;
```

4) Next, we need to convert our fasta alignments to phylip so we can run raxml. We'll use the phyluce wrappers for this (code quite liberarlly borrowed from Carl Oliveros - https://github.com/carloliveros/uce-scripts/blob/master/Species%20tree.md - thanks Carl!). Set the --alignments argument in the convert script to be where your fasta files are. You'll also want to select an appropriate outgroup for the 'phyluce_genetrees_run_raxml_genetrees.py' wrapper (if you are coming from RadSEQ/other next-gen methods, then you'll need to have phyluce installed... if you are coming from a UCE background you probably have it installed already!):

```
cd ..
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignments fasta --output phylip/ --input-format fasta --output-format phylip
```

5) We then use the phyluce wrappers to run raxml for each of our loci. Make sure to change the outgroup species in the code below to whatever it should be for your dataset 
```
python /public/uce/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_genetrees.py --input phylip --output phase_genetrees --outgroup change_this_to_your_outgroup_sp --cores 6 --quiet 2>&1 | tee kaloula_log/raxml_genetrees.txt
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



