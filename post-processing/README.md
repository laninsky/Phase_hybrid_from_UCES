#Post-processing#
Now we have our phased alleles for the hybrid tacked on to the fasta with the rest of our samples. In ONE_uce-xxx.fasta we have the first allele and the rest of our samples. In TWO_uce-xxx.fasta we have the second allele and the rest of our samples.

1) The first step is to re-align our hybrid alleles and the rest of the sequences, as indels were stripped out of the hybrid file so that gatk could eat the sequence for phasing. I am going to use MAFFT for this:
```
for i in *.fasta; do mv $i temp; mafft temp > $i; done;
```

2) Next, we need to convert our fasta alignments to phylip so we can run raxml. We'll use the phyluce wrappers for this (code quite liberarlly borrowed from Carl Oliveros - https://github.com/carloliveros/uce-scripts/blob/master/Species%20tree.md - thanks Carl!). Set the --alignments argument in the convert script to be where your fasta files are. You'll also want to select an appropriate outgroup for the 'phyluce_genetrees_run_raxml_genetrees.py' wrapper:
```
python /public/uce/phyluce/bin/align/convert_one_align_to_another.py --alignments uces --output phylip/ --input-format fasta --output-format phylip
python /public/uce/phyluce/bin/genetrees/phyluce_genetrees_run_raxml_genetrees.py --input phylip --output phase_genetrees --outgroup k_med_d --cores 6 --quiet 2>&1 | tee kaloula_log/raxml_genetrees.txt
```
