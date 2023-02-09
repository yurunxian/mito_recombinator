# Introduction

Plant mitogenomes (mitochondrial genomes) harbor variable repetitive content. These repetitive elements can sometimes mediate DNA recombination, which results in low-frequency alternative conformation of mtDNA. In the past, researchers usually used molecular experiments (i.e. PCR) to verify the existence of alternative conformations. Nowadays, application of high-throughput sequencing (especially NGS) greatly facilitates the detection of repeat-mediated recombination. By constructing alternative conformations and mapping high-throughput reads to them, one can easily confirm the recombination if the (paired) reads span the repeat and its putative flanking sequences (Fig.1).<br /><br />
![image](https://user-images.githubusercontent.com/48025559/217754280-b7a2258c-e800-424e-b443-3a0faedcc5fd.png)
Fig.1 Detection of repeat-mediated recombination with NGS reads<br /><br />

Here, I provide a pipeline to construt alternative conformations and calculate recombination frequency for each repeat.

# Pipeline

Step 1. <br />
Identify dispersed repeats for your mitogenomes with ROUSFinder.py (https://doi.org/10.25387/g3.7425680). This useful script is written in Python2 and will call BLASTN to find repeats. It is worth mentioning that its output will NOT include the contig name for each repeat. Users can manually modified the [Line:239] to make it output these information in <*_rep_table.txt> file, which is needed in the next step.

Step 2. 
