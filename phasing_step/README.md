# Phasing_step
Have you got a hybrid in your RADseq/UCE dataset, and you want to confirm its parent species? Then boy, do I have the pipeline for you! This first step ("Phase_hybrid_from_data") is pulling out the phased data for our hybrid. A subsequent step will deal with figuring out which species each allele is most closely related to.

It is assumed you have already run your data through the established RADseq/UCE initial processing pipelines (e.g. STACKS, pyRAD, https://github.com/carloliveros/uce-scripts, https://github.com/faircloth-lab/phyluce etc). For this program to work, you have access to a folder full of fasta alignments (full alignments, not just SNPs) for the loci you want to use (with any missing samples padded out with Ns or ?s), you know the samplename of your hybrid, and you've got the phasing_shell.sh script, extract_hybrid.R and onelining.R in your directory along with your phase_settings file (see below). You've got your cleaned reads in a fastq.gz file (with separate forward and reverse files if you have paired end sequencing), and you know the path to these reads.

You'll also need to have installed bwa, samtools, R and Java, and added these to your path. You'll also need to install GenomeAnalysisTK.jar (GATK) and picard.jar (picard), but we'll actually need the full pathway to these jars in the phase_settings folder below. 

###phasing_shell.sh###
The shell script is using bwa, gatk, samtools and R to pull out the hybrid sample (R), do a reference-guided assembly (bwa, samtools) on your cleaned *.fastq.gz reads from your hybrid, and then calling variants/phasing these (gatk), before using the "new reference" to do the process again to get the other alleles for your hybrid.

To run this yourself you will need a file with the input settings named phase_settings in the folder with your fasta sequences. In this file, on each separate line in this order you will need:

Line 1: the absolute pathway to GenomeAnalysisTK.jar e.g. /home/a499a400/bin/GenomeAnalysisTK.jar

Line 2: the absolute pathway to picard.jar e.g. /home/a499a400/bin/picard/dist/picard.jar

Line 3: the total number of taxa in your fasta files e.g. 24

Line 4: your hybrid sample name as it appears in your fasta files e.g. k_pixme

Line 5: either paired or single depending on your sequencing method e.g. paired

Line 6: If you had paired on the previous line, here you need the pathway to your cleaned forward reads. If you had single, then just the general pathway to your cleaned reads e.g. /home/a499a400/Kaloula/cleaned-reads/kaloula_pictameridionalishybrid_rmb586/split-adapter-quality-trimmed/kaloula_pictameridionalishybrid_rmb586-READ1.fastq.gz

Line 7: If you wrote paired on Line 5, then here you need the pathway to your cleaned reverse reads e.g. /home/a499a400/Kaloula/cleaned-reads/kaloula_pictameridionalishybrid_rmb586/split-adapter-quality-trimmed/kaloula_pictameridionalishybrid_rmb586-READ2.fastq.gz

See the example phase_settings file under 'examples'
