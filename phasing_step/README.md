# Phasing_step
Have you got a hybrid in your RADseq/UCE dataset, and you want to confirm its parent species? Then boy, do I have the pipeline for you! This first step ("Phase_hybrid_from_data") is pulling out the phased data for our hybrid. A subsequent step will deal with figuring out which species each allele is most closely related to.

#Steps you need to do before running the script#

It is assumed you have already run your data through the established RADseq/UCE initial processing pipelines (e.g. STACKS, pyRAD, https://github.com/carloliveros/uce-scripts, https://github.com/faircloth-lab/phyluce etc). For this program to work, you need to have a folder that contains:

-- fasta alignments (full alignments, not just SNPs. You also need a separate fasta file per locus, not one giant concatenated file) for the loci you want to use (with any missing samples padded out with Ns or ?s), 

-- phasing_shell.sh

-- The R-scripts: extract_hybrid.R, onelining.R and allelifying.R scripts

-- Your phasing_settings file (see below).

You'll also need to have installed bwa, samtools, R and Java, and added these to your path. You'll also need to install GenomeAnalysisTK.jar (GATK) and picard.jar (picard), but we'll actually need the full pathway to these jars in the phasing_settings folder below. 

#Coming from pyRAD#

Your .loci file can be turned in a folder of fasta alignments using the scripts at:
https://github.com/laninsky/Phase_hybrid_from_next_gen/tree/master/helper_scripts

#Getting your phasing_settings file together#

The shell script is using bwa, gatk, samtools and R to pull out the hybrid sample (R), do a reference-guided assembly (bwa, samtools) on your cleaned *.fastq.gz reads from your hybrid, and then calling variants/phasing these (gatk), before using the "new reference" to do the process again to get the other alleles for your hybrid.

To run this yourself you will need a file with the input settings named phasing_settings in the folder with your fasta sequences and the scripts. In this file, on each separate line in this order you will need:

Line 1: the absolute pathway to GenomeAnalysisTK.jar e.g. /home/a499a400/bin/GenomeAnalysisTK.jar

Line 2: the absolute pathway to picard.jar e.g. /home/a499a400/bin/picard/dist/picard.jar

Line 3: the total number of taxa in your fasta files e.g. 24

Line 4: your hybrid sample name as it appears in your fasta files e.g. k_pixme [only one "hybrid" can be processed at a time. If you are doing non-hybrids as a control, you'll also put their name here]

Line 5: either paired or single depending on your sequencing method e.g. paired

Line 6: If you had paired on the previous line, here you need the pathway to your cleaned forward reads. If you had single, then just the general pathway to your cleaned reads e.g. /home/a499a400/Kaloula/cleaned-reads/kaloula_pictameridionalishybrid_rmb586/split-adapter-quality-trimmed/kaloula_pictameridionalishybrid_rmb586-READ1.fastq.gz

Line 7: If you wrote paired on Line 5, then here you need the pathway to your cleaned reverse reads e.g. /home/a499a400/Kaloula/cleaned-reads/kaloula_pictameridionalishybrid_rmb586/split-adapter-quality-trimmed/kaloula_pictameridionalishybrid_rmb586-READ2.fastq.gz

See the example phasing_settings file under 'examples'

#To run the script#

bash phasing_shell.sh

#Is the script running successfully?#

The first few steps of phasing_shell.sh are bash/R so it might not look like it is doing much using top/htop. To confirm it is actually running, check your directory: the hybrid_ref.fa file should be growing in size as the scripts syphon off the hybrid into it. After these first steps, you should see bwa/java running through top/htop.
