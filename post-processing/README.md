#Post-processing#
Now we have our phased alleles for all our samples in separate fasta files per locus.

1) The first step is to copy all of these files to a new folder, just in case we stuff up at any of the downstream steps, so that we don't have to repeat all of the phasing steps
```
mkdir working_fasta
cp combined_fasta/*.fasta working_fasta
```

2) Next, we need to convert our fasta alignments to phylip so we can run raxml. We'll use the phyluce wrappers for this (code quite liberarlly borrowed from Carl Oliveros - https://github.com/carloliveros/uce-scripts/blob/master/Species%20tree.md - thanks Carl! If you are coming from RadSEQ/other next-gen methods, then you'll need to have python/phyluce installed... if you are coming from a UCE background you probably have it installed already!). Tweak the pathway to your convert_one_align_to_another.py file within the phyluce installation. If your sample names are longer than 10 characters you will need to pass the arguments from --shorten-names onwards. An example of the short.conf file can be found at https://github.com/laninsky/UCE_processing_steps#species-tree-stuff. Make sure you list each sample twice (once for each of its "allele names") e.g.
```
[taxa]
kaloula_baleata_jam3573_1:ba_j3573_1
kaloula_baleata_jam3573_2:ba_j3573_2
```
Run this code by:
```
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignments working_fasta --output phylip/ --input-format fasta --output-format phylip --shorten-names --name-conf  short.conf
```

5a) We can then use the phyluce wrappers to run raxml for each of our loci IF WE HAVE A COMPLETE DATASET (i.e. no missing samples for any loci - if this is you, see step 5b). You'll want to select an appropriate outgroup e.g. ba_j3573_1, tweak the number of cores if 6 is too few/too many, and edit the pathway to phyluce_genetrees_run_raxml_genetrees.py. You also need raxml to be installed and in your path.
```
python /public/uce/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_genetrees.py --input phylip --output phase_genetrees --outgroup change_this_to_your_outgroup_sp --cores 6 --quiet 2>&1 | tee log/raxml_genetrees.txt

```

5b) If we have missing samples for some of our loci, we can't use the wrapper because it assumes the outgroup we specify is present in every locus. So, we have to do this a little more manually...make sure raxml is installed and in your path. If your raxml is called something different to 'raxmlHPC-SSE3', modify this in the code below.
```
mkdir phase_genetrees
screen
```
After 'screen', you might need to cd to the working directory (the one above your phylip directory), then...
```
wd=`pwd`
cd phylip
unset i

for i in `ls *.phylip`;
do mkdir ../phase_genetrees/$i
toraxml="raxmlHPC-SSE3 -m GTRGAMMA -n best -s $wd/phylip/$i -p $RANDOM -w $wd/phase_genetrees/$i --no-bfgs"
$toraxml;
done;
```
To detach from this screen: Ctrl+A, Ctr+D. 

Note: if you get the error: ```: illegal option -- - ``` then you might need to modify the following line in the above code:
```toraxml="raxmlHPC-SSE3 -m GTRGAMMA -n best -s $wd/phylip/$i -p $RANDOM -w $wd/phase_genetrees/$i --no-bfgs"```
to specify the number of threads that are needed by -T

6) We then need to navigate to the folder with all our genetrees in it. We create a file similar to 'all-best-trees.tre', except with UCE names tied to each tree:
```
cd phase_genetrees
unset i
touch ubertree.tre
for i in `ls --color=never`; do if [ -d $i ]; then printname=`cat $i/RAxML_bestTree.best`; echo $i $printname >> ubertree.tre; fi; done 
```


#### UP TO HERE MODIFYING
3) You might be in the situation where you have multiple hybrids/samples of uncertain taxonomic status in your dataset. These are going to muddle things up a little, because if putative hybrids you have phased map back to one of these 'uncertains', that doesn't give you a lot of information on the putative parental species of the hybrid you are interested in. So... at this point we can run remove_uncertains.R before we move on to the raxml stage. To do this, you'll need to set up the file 'species_assignments'. You might as well do this now, anyway, because you are going to need it later on. Also even if you aren't planning on removing any uncertains, run the script anyway, because it helps reformat the file after MAFFT wraps the lines in it, and removes any samples completely comprised of missing data. This file has the sample names (following "samplenames.txt" from https://github.com/laninsky/reference_aligning_to_established_loci) in the first column (tab delimited), and the species/taxon designations you'd like to use in the second column. Make sure to put 'hybrid' for your putative hybrid, and 'uncertain' for any samples you think should be removed from your dataset. Also make sure to remake this file if you are running this pipeline multiple different times for different putative hybrids (or non-hybrids to test as controls). Example of 'species_assignments':
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



7) Now we are going to use some R code (hybrid_relationships.R) to pull out the clades with just our hybrid in it. This assumes you've already installed the library package stringr, data.table and plyr e.g. install.packages("stringr"). Upload hybrid_relationships.R to the folder with ubertree.tre (the 'phase_genetrees' folder), and make sure a copy of species_assignments is also in this folder (following the code below should copy species_assignments over).

```
cp ../fasta/species_assignments ../phase_genetrees
Rscript hybrid_relationships.R
```

The program will spit out a summary of relationships across your alleles for the loci. All going well, there should now be a file in your working directory called "allele_combinations.txt". You can examine it and see what combinations of closest relatives for each allele were found across the gene-trees. "NSSS" is an abbreviation for "no single sister species" e.g. in this case the hybrid allele was sister to a clade containing multiple taxa.

