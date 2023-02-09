# Introduction

Plant mitogenomes (mitochondrial genomes) harbor variable repetitive content. These repetitive elements can sometimes mediate DNA recombination, which results in low-frequency alternative conformation of mtDNA. In the past, researchers usually used molecular experiments (e.g. PCR) to verify the existence of alternative conformations. Nowadays, application of high-throughput sequencing (especially NGS) greatly facilitates the detection of repeat-mediated recombination. By constructing alternative conformations and mapping high-throughput reads to them, one can easily confirm the recombination if the (paired) reads span the repeat and its putative flanking sequences (Fig.1).<br /><br />
![image](https://user-images.githubusercontent.com/48025559/217754280-b7a2258c-e800-424e-b443-3a0faedcc5fd.png)
Fig.1 Detection of repeat-mediated recombination with NGS reads<br /><br />

Here, I provide a pipeline to construct alternative conformations and calculate recombination frequency for plant mitogenomes.

# Pipeline

Step 1. <br />
Identify dispersed repeats for your mitogenomes with ROUSFinder.py (https://doi.org/10.25387/g3.7425680). This useful script is written in Python2 and will call BLASTN to find repeats. It is worth mentioning that its output will NOT include the contig name for each repeat. Users can manually modified the [Line:239] to make it output these information in <*_rep_table.txt> file, which is needed in the next step.

Step 2. <br />
Use **recombination.pl** to construct four conformations for each repeat (two genomic and two alternative). Each conformation contains the repetitive sequence itself and 300 bp flanking sequences. <br />
`perl recombination.pl <*_rep_table.txt> <mitogenome.fasta> <conformation.txt> <conformation.statistics.txt>`<br />
This script only process repeat pair (two-copy sequence) that meets the following standards:
> (1) There is no overlap and a >50 bp distance between copies.<br />
> (2) At least one copy has no overlap with another larger repetitive sequence.<br />
> (3) To avoid the influence of nested tandem repeats (especially for large repeat), for each alternative conformation, the core repetitive sequence and its [half of the insert size for paired reads (350 bp in this case)] bp flanking sequences should not have an identical counterpart in the other alternative conformation.

Step 3. <br />
Map your mitochondrial reads (of course paired) to <conformation.txt> and get the sam file.<br />
`bowtie2-build <conformation.txt> conformation`<br />
`bowtie2 -x conformation -1 <left.fq> -2 <right.fq> --end-to-end --no-discordant --no-mixed --very-fast -p 20 --no-unal -S <conformation.sam>`<br />
Extract the key information from your sam file.<br />
`cut -f1,3,4,6,9 <conformation.sam> <conformation.txt>`<br />
Alternatively, one may want you use my **mapping.pl**<br />.
`perl mapping.pl <conformation.txt> <conformation.statistics.txt> <left.fq> <right.fq> <conformation.sam>`

Step 4. <br />
Use **recombination_ratio.pl** to calculate recombination frequency for each repeat.<br />
`perl recombination_ratio.pl <conformation.statistics.txt> <conformation.txt> <conformation.txt> <conformation.statistics.done.txt>`<br />
Then you can find your results in <conformation.statistics.done.txt>.

![}HHX_NA%02Q7}HDXW@LDSNS](https://user-images.githubusercontent.com/48025559/217771027-1ef4b9d9-7d24-4715-8001-d8662d136c5e.png)

