#!/bin/sh

sample=$1

for lib in `seq 1 4`
do bowtie2 -x ~/Genomes/Bowtie2_indexes/Branchiostoma_lanceolatum/B_lanceolatum -U Reads_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_in_$sample'_'$lib'.fastq' --no-unal --quiet -S Genome_mapping_of_non_ncRNA_non_pre-miRNA_matching_reads_$sample'_'$lib'.sam'
done
for lib in `seq 1 4`
do sed '1,/^@PG\t/ !d' Genome_mapping_of_non_ncRNA_non_pre-miRNA_matching_reads_$sample'_'$lib'.sam' > 18_to_30-mers_Genome_mapping_of_non_ncRNA_non_pre-miRNA_matching_reads_$sample'_'$lib'.sam'
   sed '1,/^@PG\t/ d' Genome_mapping_of_non_ncRNA_non_pre-miRNA_matching_reads_$sample'_'$lib'.sam' | awk 'length($10)>=18 && length($10)<=30 {print}' >> 18_to_30-mers_Genome_mapping_of_non_ncRNA_non_pre-miRNA_matching_reads_$sample'_'$lib'.sam'
done
for lib in `seq 1 4`
do samtools view -Sb 18_to_30-mers_Genome_mapping_of_non_ncRNA_non_pre-miRNA_matching_reads_$sample'_'$lib'.sam' | bedtools intersect -bed -wa -a stdin -b RdRP_template_loci.bed > genomic_mapping_of_$sample'_'$lib'_on_RdRP_template_loci.bed'
done

