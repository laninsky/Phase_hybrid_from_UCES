# Phase_hybrid_from_next-gen_data
Have you got a hybrid in your next-gen dataset (currently developed for UCE and RadSEQ loci), and you want to confirm its parent species? Then boy, do I have the pipeline for you! You can either manually copy the R and bash scripts to your local computer, or use git clone.

There are two sets of scripts to help with this (their relative contribution is outlined in the image below):

The first step is in the folder phasing_step. There is a readme in that folder to help with those.

The second step is in the folder post-processing. There is also a readme in there to help you too!

If you are coming from pyRAD, you can use the scripts in helper_scripts to get your fasta files together from the *.loci file.

You'll probably want to run the pipeline more than once: for each of your hybrids, and for a few non-hybrids to act as a comparison. If you are interested in just going ahead and phasing your full dataset, see https://github.com/laninsky/phase_everyone

![hybrid_identification_pipeline](https://cloud.githubusercontent.com/assets/8808649/8860428/4e512bf6-3149-11e5-9ee4-7ea321904d3b.png)

#Suggested citation
This code was first published in: TBD

If you could cite the pub, and the progam as below, that would be lovely:

Alexander, A. 2015. Phase_hybrid_from_next_gen v1.0.0. Available from https://github.com/laninsky/Phase_hybrid_from_next_gen

#Version history
v1.0.0: ready to rock 'n' roll
