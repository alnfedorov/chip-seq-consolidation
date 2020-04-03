## What is the pain the pain to address?
How to reliably merge several biological replicates? Especially if they don't have a lot of reads.
* IDR - only two replicates in practice, bad for broad histon marks
* Pool - doesn't account for biases in the specific platforms experiments. Might be good in some cases.
* Peak-calling and merge later. Good option, but how reliable is it? How to merge?

For ULI CHip-seq it is common to have few reads for the peak-source histone marks. Hard to reliably peak-call ~10m reads. 

CHip-seq has it's own problems - antibody quality, sequencing and alignment problems. Errors on each step introduces 
unnecessary biases that affect peak-calling algorithms.

What do I want? An option to reliably merge several biological replicates. Even if there is problems during peak-calling 
or any other pipeline steps. How to do this? Peaks MUST have some sort of biological similarity inside/around them in order to 
distinguishable from other genomic locations. Hence we can relay on the ability of neural networks to learn in presence 
of the noisy labels to filter out false-positive and probably false-negative regions.
## QC options
### Tuning sandbox
##### Simulations:  
Simulate chip-seq experiment from scratch.  
* Requires high-quality gt regions(not available for histone's chip-seq).  
* Sophisticated and perhaps unreliable tooling.  
##### Real data:
Use ENCODE data.  
* Low-quality and noisy.  
* There is no gt regions. Hence, it is hard to measure the quality of the reconstruction.  
##### Hybrid approach
Use deeply sequenced epigenomes as a golden standart and apply several peak-calling tools to obtain reliable set of 
peaks. Use them as a ground true for simulations.  
* More or less reliable approach.  
* Use simulations with random set of peaks as a check for adequacy in simulations.  
### Metrics
1. IOU between replicates and consensus.
2. Precision and recall for sandbox replicas.
3. Fraction of reads in peaks?
4. Signal to noise ratio?
# TODO
1. Tooling to download/validate and concatenate fastq data from ENCODE
2. Tooling to align(bowtie2), call-peaks(macs2(narrow-broad), epic2-sicer(narrow-broad), 
                                         peakranger bcp(broad), peakranger ranger(narrow),
                                         PePr(narrow-broad), MUSIC(narrow-broad))
3. Tooling for each peak-calling util

What steps do I need to cover?
2. Ability to parse arbitrary encode experiment(later)
3. Utils to run peak-calling
4. Utils to run bowtie2 aligning 
5. End-to-end pipeline script

TODO:
1. Make peak-calling for the sake of humanity.


