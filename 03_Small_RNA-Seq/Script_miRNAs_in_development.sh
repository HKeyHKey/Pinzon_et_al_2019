#!/bin/sh

TOLERANCE_3p=3 # accepts up to $TOLERANCE_3p trimmed nucleotides on the 3´ end (and will accept any length of templated 3´ extension)
ABUNDANCE_CUTOFF_FOR_REPORT=10 # minimal number of ppm (in at least 1 developmental stage) for the (5´ arm or 3´ arm) miRNA to be reported in the table

echo "Hairpin_ID Sample Reads_for_5p_RNA Reads_for_3p_RNA" > miRNA_abundance_in_development.dat

echo "\\\begin{longtable}{|c|c|c|}" > table_miRNAs_development.tex
echo "\\hline" >> table_miRNAs_development.tex
echo "Pre-miRNA&miRNA sequences&Abundance profile in development\\\\\\" >> table_miRNAs_development.tex
echo "\\hline" >> table_miRNAs_development.tex
echo "\\endfirsthead" >> table_miRNAs_development.tex
echo "\\hline" >> table_miRNAs_development.tex
echo "Pre-miRNA&miRNA sequences&Abundance profile in development\\\\\\" >> table_miRNAs_development.tex
echo "\\hline" >> table_miRNAs_development.tex
echo "\\endhead" >> table_miRNAs_development.tex
echo "\\hline \\multicolumn{3}{|r|}{\\scriptsize{(continued on next page)}\\\\normalsize}\\\\\\ \\hline" >> table_miRNAs_development.tex
echo "\\endfoot" >> table_miRNAs_development.tex
echo "\\hline" >> table_miRNAs_development.tex
echo "\\endlastfoot" >> table_miRNAs_development.tex
echo "\\hline" >> table_miRNAs_development.tex

for hairpin in `grep '^>' Convincing_pre-miRNA_hairpins.fa | sed -e 's|^> *||' -e 's| .*||'`
do loop_size=`grep '^'$hairpin' ' Hairpin_structure.dat | awk '{print $3}' | grep -o '(\.*)' | wc -m`
   loop_size=`echo $loop_size"-3" | bc` # to subtract the apical base pair and the carriage return
   last_pair=`grep '^'$hairpin' ' Hairpin_structure.dat | awk '{print $3}' | sed 's|(\.\.*).*||' | wc -m`
   middle_loop=`echo "scale=1;"$last_pair"+"$loop_size"/2" | bc`
   for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male
   do sed '1,/^@PG\t/ d' Convincing_hairpin-matching_reads_$sample'_1.sam' | awk '$2==0 && $3=='$hairpin' {print $4,$4+length($10)-1}' > tmp_hits_$hairpin'_'$sample
      cat tmp_hits_$hairpin'_'$sample
   done | sort | uniq -c > tmp_total_hits_$hairpin
   major_5p=`awk '$3<'$middle_loop' {print}' tmp_total_hits_$hairpin | sort -g | tail -1 | awk '{print $2,$3}'`
   major_3p=`awk '$2>'$middle_loop' {print}' tmp_total_hits_$hairpin | sort -g | tail -1 | awk '{print $2,$3}'`
   start_major_5p=`echo $major_5p | awk '{print $1}'`
   end_major_5p=`echo $major_5p | awk '{print $2}'`
   start_major_3p=`echo $major_3p | awk '{print $1}'`
   end_major_3p=`echo $major_3p | awk '{print $2}'`
   for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male 
   do abundance_5p=`awk '$1=='$start_major_5p' && $2>='$end_major_5p'-'$TOLERANCE_3p' {print}' tmp_hits_$hairpin'_'$sample | wc -l`
      abundance_3p=`awk '$1=='$start_major_3p' && $2>='$end_major_3p'-'$TOLERANCE_3p' {print}' tmp_hits_$hairpin'_'$sample | wc -l`
      echo $hairpin $sample $abundance_5p $abundance_3p >> miRNA_abundance_in_development.dat
   done

   hairpin_seq=`grep -A 1 '^> *'$hairpin' ' Convincing_pre-miRNA_hairpins.fa | tail -1`
   major_5p=`echo $major_5p | sed 's| |-|'`
   major_3p=`echo $major_3p | sed 's| |-|'`
   miRNA_5p_seq=`echo $hairpin_seq | cut -c $major_5p | tr T U`
   miRNA_3p_seq=`echo $hairpin_seq | cut -c $major_3p | tr T U`

### Below: tests whether the 5´ arm and 3´ arm miRNAs are abundant enough to be worth reporting:
   ab5p=0
   ab3p=0
   for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male
   do raw_5p=`awk '$1=='$hairpin' && $2=="'$sample'" {print $3}' miRNA_abundance_in_development.dat`
      raw_3p=`awk '$1=='$hairpin' && $2=="'$sample'" {print $4}' miRNA_abundance_in_development.dat`
      depth=`awk '$1=="'$sample'" && $2==1 {print $5-$6}' Statistics.dat`
      normalized_5p=`echo $raw_5p"*1000000/"$depth | bc`
      normalized_3p=`echo $raw_3p"*1000000/"$depth | bc`
      if test $normalized_5p -gt $ABUNDANCE_CUTOFF_FOR_REPORT
      then ab5p=1
      fi
      if test $normalized_3p -gt $ABUNDANCE_CUTOFF_FOR_REPORT
      then ab3p=1
      fi
   done

#   echo "raw_5p="$raw_5p" raw_3p="$raw_3p" ab5p="$ab5p" ab3p="$ab3p

   if test $ab5p -eq 0
   then miRNA_5p_seq='(low abundance)'
   fi
   if test $ab3p -eq 0
   then miRNA_3p_seq='(low abundance)'
   fi  
### Above: tests whether the 5´ arm and 3´ arm miRNAs are abundant enough to be worth reporting

   annot=`grep -B 1 -w -m 1 $hairpin_seq Orthologs_to_known_hairpins.fa | grep '^>' | sed 's|^> *||'`
   annot1=`echo $annot | sed 's| (.*||'`
   annot2=`echo $annot | sed 's|.* (|(|'`
   echo "\\\begin{minipage}[b]{2.1cm} "$annot1"\\\\\\"$annot2"\\end{minipage}" >> table_miRNAs_development.tex
   echo "&" >> table_miRNAs_development.tex
   echo "\\\begin{minipage}[b]{4.6cm} 5\\'{}arm:\\\\\\\\\\\begin{complem}"$miRNA_5p_seq"\\end{complem}" >> table_miRNAs_development.tex
   echo "" >> table_miRNAs_development.tex
   echo "\\\bigskip" >> table_miRNAs_development.tex
   echo "" >> table_miRNAs_development.tex
   echo "3\\'{}arm:\\\\\\\\\\\begin{complem}"$miRNA_3p_seq"\\end{complem}" >> table_miRNAs_development.tex
   echo "\\end{minipage}" >> table_miRNAs_development.tex
   echo "&" >> table_miRNAs_development.tex
   echo "\\includegraphics[width=10cm]{/home/herve/Amphioxus/miRNA_search/Developmental_profile_miRNAs_from_hairpin_"$hairpin".pdf}\\\\\\" >> table_miRNAs_development.tex
   echo "\\hline" >> table_miRNAs_development.tex
done
echo "\\end{longtable}" >> table_miRNAs_development.tex
R CMD BATCH R_commands_plot_miRNAs_in_dvpt
