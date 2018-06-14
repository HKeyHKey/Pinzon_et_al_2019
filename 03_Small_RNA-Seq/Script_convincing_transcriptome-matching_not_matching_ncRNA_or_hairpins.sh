#!/bin/sh

sample=$1
for lib in `seq 1 4`
do bowtie2 -x Convincing_ORFs -U Reads_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_in_$sample'_'$lib'.fastq' --no-unal --quiet -S Convincing_ORF_transcriptome-matching_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_$sample'_'$lib'.sam'
done
./Chef4.sh $sample
