#Post-processing#
Now we have our phased alleles for all our samples in separate fasta files per locus.

1) The first step is to copy all of these files to a new folder, just in case we stuff up at any of the downstream steps, so that we don't have to repeat all of the phasing steps
```
mkdir working_fasta
cp combined_fasta/*.fasta working_fasta
```

2) Next, we need to convert our fasta alignments to phylip so we can run raxml. We'll use the phyluce wrappers for this (code quite liberally borrowed from Carl Oliveros - https://github.com/carloliveros/uce-scripts/blob/master/Species%20tree.md - thanks Carl! If you are coming from RadSEQ/other next-gen methods, then you'll need to have python/phyluce installed... if you are coming from a UCE background you probably have it installed already!). Tweak the pathway to your convert_one_align_to_another.py file within the phyluce installation. If your sample names are longer than 10 characters you will need to pass the arguments from --shorten-names onwards. An example of the short.conf file can be found at https://github.com/laninsky/UCE_processing_steps#species-tree-stuff. Make sure you list each sample twice (once for each of its "allele names") e.g.
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

7) Now we are going to use some R code (hybrid_relationships.R) to pull out the relationships of the alleles for each of our samples. This assumes you've already installed the library package stringr, data.table and plyr e.g. install.packages("stringr"). Upload hybrid_relationships.R to the folder with ubertree.tre (the 'phase_genetrees' folder). You'll also need a file in the same folder which has the allele assignment suffix for allele 1 in the first row, and the assignment suffix for allele 2 in the second row called "suffix.txt" e.g.
```
_1
_2
```
Then invoke the script by:
```
Rscript hybrid_relationships.R
```

The program will spit out a summary of relationships across the alleles for all loci for each sample. All going well, there should now be a file in your working directory called "allele_combinations.txt". You can examine it and see what combinations of closest relatives for each allele were found for each sample at specific loci. "allele_combinations_sum.txt" has the same information but summarized for each sample. Relationships are not recorded for alleles which are sister to a clade of other alleles rather than to a single other allele.
