#!/bin/sh

sample=$1

for lib in `seq 1 4`
do sed '1,/^@PG\t/ !d' Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam' > 18_to_30-mers_Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam'
   sed '1,/^@PG\t/ d' Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam' | awk 'length($10)>=18 && length($10)<=30 {print}' >> 18_to_30-mers_Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam'
done

for lib in `seq 1 4`
do samtools view -Sb 18_to_30-mers_Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam' | bedtools intersect -bed -wa -a stdin -b RdRP_template_transcripts.bed > extragenomic_mapping_of_$sample'_'$lib'_on_RdRP_template_transcripts.bed'
done

