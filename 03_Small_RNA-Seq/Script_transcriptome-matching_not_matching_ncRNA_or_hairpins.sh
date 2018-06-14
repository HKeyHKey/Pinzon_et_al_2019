#!/bin/sh

sample=$1
for lib in `seq 1 4`
do cat /mnt/data/home/herve.seitz/Amphioxus/Small_RNA_mapping/All_reads_in_$sample'_'$lib'.txt' /mnt/data/home/herve.seitz/Amphioxus/Small_RNA_mapping/IDs_matching_abundant_ncRNAs_$sample'_'$lib'.txt' /mnt/data/home/herve.seitz/Amphioxus/miRNA_search/IDs_matching_hairpins_$sample'_'$lib'.txt' | sort | uniq -c | grep '^ *1 ' | awk '{print $2}' > Reads_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_in_$sample'_'$lib'.txt'
   ./Module_extract_from_fastq.pl /mnt/data/home/herve.seitz/Natalia/projet_Amphioxus/$sample'_'$lib'_Trimmed'* Reads_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_in_$sample'_'$lib'.txt' > Reads_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_in_$sample'_'$lib'.fastq'
   bowtie2 -x /mnt/data/home/herve.seitz/Genomes/Bowtie2_indexes/Branchiostoma_lanceolatum/B_lanceolatum_transcriptome -U Reads_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_in_$sample'_'$lib'.fastq' --no-unal --quiet -S Transcriptome-matching_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_$sample'_'$lib'.sam'
done
