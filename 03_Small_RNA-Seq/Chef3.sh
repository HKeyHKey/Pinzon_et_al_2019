#!/bin/sh

sample=$1

for lib in `seq 1 4`
do sed '1,/^@PG\t/ d' Transcriptome-matching_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_$sample'_'$lib'.sam' | awk '$2==0 {print ">"$1"\n"$10}' > sense_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'$lib'.fa'
   sed '1,/^@PG\t/ d' Transcriptome-matching_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_$sample'_'$lib'.sam' | awk '$2==16 {print ">"$1}' > tmp1_$sample
   sed '1,/^@PG\t/ d' Transcriptome-matching_not_matching_abundant_ncRNAs_nor_pre-miRNA_hairpins_$sample'_'$lib'.sam' | awk '$2==16 {print $10}' | rev | tr ACGT TGCA > tmp2_$sample
   paste -d "\n" tmp1_$sample tmp2_$sample > antisense_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'$lib'.fa'
   sed -i '/^>/ !s|T|U|g' sense_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'$lib'.fa'
   sed -i '/^>/ !s|T|U|g' antisense_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'$lib'.fa'
   echo "Size Sense_matching Antisense_matching" > Size_distribution_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'$lib'.dat'
   for size in `seq 18 30`
   do display=$size
      for ori in sense antisense
      do grep -B 1 '^[ACGU]\{'$size'\}$' $ori'_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_'$sample'_'$lib'.fa' | grep -v '\-\-' > $size'-mers_'$ori'_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_'$sample'_'$lib'.fa'
         nb_seq=`grep -c '>' $size'-mers_'$ori'_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_'$sample'_'$lib'.fa'`
         display=`echo $display $nb_seq`
         if test -s $size'-mers_'$ori'_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_'$sample'_'$lib'.fa'
         then /mnt/data/home/herve.seitz/Annotation_pipeline/weblogo/seqlogo -f $size'-mers_'$ori'_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_'$sample'_'$lib'.fa' -F EPS -o Logo_$size'-mers_'$ori'_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_'$sample'_'$lib -abcMnY
         fi
      done
      echo $display >> Size_distribution_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'$lib'.dat'
  done
done
tar -cf Figure_3-style_$sample'.tar' Size_distribution_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'* Logo_*_transcriptome-matching_reads_not_matching_abundant_ncRNAs_or_pre-miRNA_hairpins_$sample'_'*;bzip2 Figure_3-style_$sample'.tar'
