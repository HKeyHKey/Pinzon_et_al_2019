#!/bin/sh

sample=$1

for lib in `seq 1 4`
do for species in ATCV1
   do sed '1,/^@PG\t/ d' $species'_mapping_'$sample'_'$lib'.sam' | awk '{print ">"$1"\n"$10}' | sed '/^>/ !s|T|U|g' > Extragenomic_extratranscriptomic_reads_matching_$species'_in_'$sample'_'$lib'.fa'
      echo "Size Number_of_reads" > Size_distribution_$species'_reads_'$sample'_'$lib'.dat'
      for size in `seq 18 30`
      do display=$size
         grep -B 1 '^[ACGU]\{'$size'\}$' Extragenomic_extratranscriptomic_reads_matching_$species'_in_'$sample'_'$lib'.fa' | grep -v '\-\-' > $size'-mers_Extragenomic_extratranscriptomic_reads_matching_'$species'_in_'$sample'_'$lib'.fa'
         nb_seq=`grep -c '>' $size'-mers_Extragenomic_extratranscriptomic_reads_matching_'$species'_in_'$sample'_'$lib'.fa'`
         display=`echo $display $nb_seq`
         if test -s $size'-mers_Extragenomic_extratranscriptomic_reads_matching_'$species'_in_'$sample'_'$lib'.fa'
         then /mnt/data/home/herve.seitz/Annotation_pipeline/weblogo/seqlogo -f $size'-mers_Extragenomic_extratranscriptomic_reads_matching_'$species'_in_'$sample'_'$lib'.fa' -F EPS -o Logo_$size'-mers_Extragenomic_extratranscriptomic_reads_matching_'$species'_in_'$sample'_'$lib -abcMnY
         fi
         echo $display >> Size_distribution_$species'_reads_'$sample'_'$lib'.dat'
     done
   done
done
