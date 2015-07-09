#Post-processing#
Now we have our phased alleles for the hybrid tacked on to the fasta with the rest of our samples. In ONE_uce-xxx.fasta we have the first allele and the rest of our samples. In TWO_uce-xxx.fasta we have the second allele and the rest of our samples.

1) The first step is to  re-align our hybrid alleles and the rest of the sequences, as indels were stripped out of the hybrid file so that gatk could eat the sequence for phasing. I am going to use MAFFT for this, but you could pick your favorite aligher:
```
for i in *.fasta; do mv $i temp; mafft temp > $i; done;
rm temp
```

2) Now we are going to copy all of these files to a new folder, just in case we stuff up at any of the downstream steps, so that we don't have to repeat everything
```
mkdir fasta
cp *.fasta fasta
cd fasta
```

3) You might be in the situation where you have multiple hybrids/samples of uncertain taxonomic status in your dataset. These are going to muddle things up a little, because if the hybrid you have phased maps back to one of these 'uncertains', that doesn't give you a lot of information on the putative parental species of the hybrid you are interested in. So... at this point we can run remove_uncertains.R before we move on to the raxml stage. To do this, you'll need to set up the file 'species_assignments'. You might as well do this now, anyway, because you are going to need it later on. This file has the sample names in the first column (tab delimited), and the species/taxon designations you'd like to use in the second column. Make sure to put 'hybrid' for your putative hybrid, and 'uncertain' for any samples you think should be removed from your dataset. Also make sure to remake this file if you are running this pipeline multiple different times for different putative hybrids! Example of 'species_assignments':
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
unset i
for i in `ls *.fasta`;
do mv $i temp;
Rscript remove_uncertains.R
mv temp.fa $i;
rm -rf temp;
done;
```

Even if you aren't planning on removing any uncertains, run the script anyway, because it helps reformat the file after MAFFT wraps the lines in it.

4) Next, we need to convert our fasta alignments to phylip so we can run raxml. We'll use the phyluce wrappers for this (code quite liberarlly borrowed from Carl Oliveros - https://github.com/carloliveros/uce-scripts/blob/master/Species%20tree.md - thanks Carl! If you are coming from RadSEQ/other next-gen methods, then you'll need to have python/phyluce installed... if you are coming from a UCE background you probably have it installed already!). Tweak the pathway to your convert_one_align_to_another.py file within the phyluce installation:

```
cd ..
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignments fasta --output phylip/ --input-format fasta --output-format phylip
```

5) We then use the phyluce wrappers to run raxml for each of our loci. You'll want to select an appropriate outgroup e.g. k_bal_r, tweak the number of cores if 6 is too few/too many, and edit the pathway to phyluce_genetrees_run_raxml_genetrees.py
```
python /public/uce/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_genetrees.py --input phylip --output phase_genetrees --outgroup change_this_to_your_outgroup_sp --cores 6 --quiet 2>&1 | tee kaloula_log/raxml_genetrees.txt

```

6) We then need to navigate to the folder with all our genetrees in it. We create a file similar to 'all-best-trees.tre', except with UCE names tied to each tree:
```
cd phase_genetrees
unset i
touch ubertree.tre
for i in `ls --color=never`; do if [ -d $i ]; then printname=`cat $i/RAxML_bestTree.best`; echo $i $printname >> ubertree.tre; fi; done 
```

7) Now we are going to use some R code (hybrid_relationships.R) to pull out the clades with just our hybrid in it. This assumes you've already installed the library package stringr, data.table and plyr e.g. install.packages("stringr"). Upload hybrid_relationships.R to the folder with ubertree.tre (the 'phase_genetrees' folder), and make sure a copy of species_assignments is also in this folder (following the code below should copy species_assignments over).

```
cp ../fasta/species_assignments ../phase_genetrees
Rscript hybrid_relationships.R
```

The program will spit out a summary of relationships across your alleles for the loci. All going well, there should now be a file in your working directory called "allele_combinations.txt". You can examine it and see what combinations of closest relatives for each allele were found across the gene-trees. "NSSS" is an abbreviation for "no single sister species" e.g. in this case the hybrid allele was sister to a clade containing multiple taxa.

