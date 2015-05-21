# Phase_hybrid_from_UCES
Have you got a hybrid in your ultra-conserved element dataset, and you want to confirm its parent species? Then boy, do I have the pipeline for you!

It is assumed you have already run your data through the established UCE pipelines (e.g. https://github.com/carloliveros/uce-scripts, https://github.com/faircloth-lab/phyluce). Here, you have access to a folder full of fasta alignments for all your UCEs, and you know what you have named your hybrid, and you've uploaded the dependent Rscripts into the same directory full of your fasta alignments.

###phasing_shell.sh###
The shell script is using bwa, gatk, samtools and R to pull out the hybrid sample (R), do a reference-guided assembly (bwa, samtools) using cleaned-data you trimmed during the uce pipeline step, and then calling variants/phasing these (gatk), before using the "new reference" to do the process again to get the other allele for your hybrid.

Things you'll need to change in the shell script to run this yourself:

--the pathways to the programs

--the pathways to your reads

--any mention of 'pixme' (my hybrid species!) needs to be replaced with your own hybrid species name.

--you'll also need to change some of this stuff inside the R code that this shell calls (extract_hybrid.R and onelining.R)

###extract_hybrid.R####


###onelining.R###
