#!/bin/sh

sample=$1
contig=$2

for lib in `seq 1 4`
do sed '1,/^@PG\t/ d' Unified_contig_mapping_$sample'_'$lib'.sam' | awk '$3=="'$contig'" {print ">"$1"\n"$10}' | sed '/^>/ !s|T|U|g' > $contig'_matching_in_'$sample'_'$lib'.fa'
   echo "Size Number_of_reads" > Size_distribution_$contig'_reads_'$sample'_'$lib'.dat'
   for size in `seq 18 30`
   do display=$size
      grep -B 1 '^[ACGU]\{'$size'\}$' $contig'_matching_in_'$sample'_'$lib'.fa' | grep -v '\-\-' > $size'-mers_'$contig'_in_'$sample'_'$lib'.fa'
      nb_seq=`grep -c '>' $size'-mers_'$contig'_in_'$sample'_'$lib'.fa'`
      display=`echo $display $nb_seq`
      if test -s $size'-mers_'$contig'_in_'$sample'_'$lib'.fa'
      then /mnt/data/home/herve.seitz/Annotation_pipeline/weblogo/seqlogo -f $size'-mers_'$contig'_in_'$sample'_'$lib'.fa' -F EPS -o Logo_$size'-mers_'$contig'_in_'$sample'_'$lib -abcMnY
      fi
      echo $display >> Size_distribution_$contig'_reads_'$sample'_'$lib'.dat'
   done
done
