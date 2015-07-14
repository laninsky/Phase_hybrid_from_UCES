# Converting pyRAD *.loci file to a whole ton of fasta files#

What you'll need in your folder:

-- your *.loci file from your pyRAD output.

-- a 'species_assignments' file (see https://github.com/laninsky/Phase_hybrid_from_next_gen/blob/master/post-processing/README.md)

-- the converting_pyrad.R script

To convert the file run the following command:

for i in 'ls *.loci';
do cp $i temp
Rscript converting_pyrad.R;
done;

Then you are ready to copy the other scripts you need for the phasing_step to this folder
