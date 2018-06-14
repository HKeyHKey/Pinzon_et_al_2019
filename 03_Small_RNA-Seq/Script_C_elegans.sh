#for f in `ls Logo_*-mers_extratranscriptomic_extragenomic_reads_*.pdf`
#do ln -s $f `echo $f | sed 's|extratranscriptomic_extragenomic|extragenomic_extratranscriptomic|'`
#done
#for f in `ls Logo_*sense_mapping_transcriptomic_extragenomic_reads_*.pdf`
#do ln -s $f `echo $f | sed 's|mapping_transcriptomic_extragenomic|transcriptomic_extragenomic|'`
#done
#for f in `ls Logo_*sense_mapping_Convincing_ORF_transcript_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_*.pdf`
#do ln -s $f `echo $f | sed 's|mapping_Convincing_ORF_transcript_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins|convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads|'`
#done
#for f in `ls Logo_*sense_mapping_Transcriptome_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_*.pdf`
#do ln -s $f `echo $f | sed 's|mapping_Transcriptome_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins|transcriptomic_not_ncRNA_not_hairpin-matching_reads|'`
#done
#for f in `ls Logo_*sense_Read_IDs_hairpin_mapping_*.pdf`
#do ln -s $f `echo $f | sed 's|Read_IDs_hairpin_mapping|hairpin-matching_reads|'`
#done
#for f in `ls Logo_*-mers_genomic_reads_not_matching_abundant_cel_ncRNA_*.pdf`
#do ln -s $f `echo $f | sed 's|genomic_reads_not_matching_abundant_cel_ncRNA_|genomic_non_ncRNA-matching_reads_|'`
#done


for analysis in genomic_non_ncRNA-matching_reads hairpin-matching_reads transcriptomic_not_ncRNA_not_hairpin-matching_reads convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads extragenomic_extratranscriptomic_reads transcriptomic_extragenomic_reads
do case "$analysis" in "genomic_non_ncRNA-matching_reads") title="Genomic reads not matching abundant ncRNAs";;
		       "hairpin-matching_reads") title="pre-miRNA hairpin-matching reads";;
                       "transcriptomic_not_ncRNA_not_hairpin-matching_reads") title="Transcriptome-matching reads (excuding pre-miRNA and abundant ncRNA-matching reads)";;
                       "convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads") title="Reads matching RNAs with long ORF's";;
                       "extragenomic_extratranscriptomic_reads") title="Extragenomic and extratranscriptomic reads";;
		       "transcriptomic_extragenomic_reads") title="Extragenomic reads matching the transcriptome";;
   esac
   echo "\\subsection{"$title"}"
   for lib in GSM455391 GSM455392 GSM455393
   do case "$lib" in "GSM455391") lib_title="18-26-mers, any number of 5\'{} phosphates, replicate 1";;
                     "GSM455392") lib_title="18-26-mers, any number of 5\'{} phosphates, replicate 2";;
                     "GSM455393") lib_title="18-26-mers, any number of 5\'{} phosphates, replicate 3";;
      esac
      echo "Library "$lib" ("$lib_title"):"
      echo ""
      echo ""
      echo "\\begin{center}"
      echo "\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Size_distribution_"$analysis"_"$lib".pdf}"
      echo "\\end{center}"
      echo ""
      if test "$analysis" = "extragenomic_extratranscriptomic_reads" -o "$analysis" = "genomic_non_ncRNA-matching_reads"
      then echo "\\begin{tabular}{ccc}"
           size1=18
           while test $size1 -le 24
           do size2=`echo $size1"+1" | bc`
              size3=`echo $size2"+1" | bc`
              echo $size1"-mers:&"$size2"-mers:&"$size3"-mers:\\\\"
              echo "\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size1"-mers_"$analysis"_"$lib".pdf}&\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size2"-mers_"$analysis"_"$lib".pdf}&\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size3"-mers_"$analysis"_"$lib".pdf}\\\\"
              size1=`echo $size1"+3" | bc`
           done
           echo "\\end{tabular}"
           echo "\\newpage"
      else echo "Sense reads:"
           echo ""
           echo "\\begin{tabular}{ccc}"
           size1=18
           while test $size1 -le 24
           do size2=`echo $size1"+1" | bc`
              size3=`echo $size2"+1" | bc`
              echo $size1"-mers:&"$size2"-mers:&"$size3"-mers:\\\\"
              echo "\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size1"-mers_sense_"$analysis"_"$lib".pdf}&\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size2"-mers_sense_"$analysis"_"$lib".pdf}&\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size3"-mers_sense_"$analysis"_"$lib".pdf}\\\\"
              size1=`echo $size1"+3" | bc`
           done
           echo "\\end{tabular}"
           echo "\\newpage"
           echo "Antisense reads:"
           echo ""
           echo "\\begin{tabular}{ccc}"
           size1=18
           while test $size1 -le 24
           do size2=`echo $size1"+1" | bc`
              size3=`echo $size2"+1" | bc`
              echo $size1"-mers:&"$size2"-mers:&"$size3"-mers:\\\\"
              echo "\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size1"-mers_antisense_"$analysis"_"$lib".pdf}&\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size2"-mers_antisense_"$analysis"_"$lib".pdf}&\\includegraphics[width=6cm]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_"$size3"-mers_antisense_"$analysis"_"$lib".pdf}\\\\"
              size1=`echo $size1"+3" | bc`
           done
           echo "\\end{tabular}"
      fi             
      echo "\\newpage"
   done
done > C_elegans.tex

for analysis in genomic_non_ncRNA-matching_reads hairpin-matching_reads transcriptomic_not_ncRNA_not_hairpin-matching_reads convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads extragenomic_extratranscriptomic_reads transcriptomic_extragenomic_reads
do for size in `seq 18 26`
   do for lib in GSM455391 GSM455392 GSM455393
      do if [ ! -f Logo_"$size"-mers_"$analysis"_"$lib".pdf ]
         then sed -i 's|\\includegraphics\[width=6cm\]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_'$size'-mers_'$analysis'_'$lib'.pdf}|(no read)|g' C_elegans.tex
         fi
	 if [ ! -f Logo_"$size"-mers_sense_"$analysis"_"$lib".pdf ]
	 then sed -i 's|\\includegraphics\[width=6cm\]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_'$size'-mers_sense_'$analysis'_'$lib'.pdf}|(no read)|g' C_elegans.tex
         fi
         if [ ! -f Logo_"$size"-mers_antisense_"$analysis"_"$lib".pdf ]
         then sed -i 's|\\includegraphics\[width=6cm\]{/home/herve/Amphioxus/Size_dist_and_logo_September2017/Logo_'$size'-mers_antisense_'$analysis'_'$lib'.pdf}|(no read)|g' C_elegans.tex
         fi
      done
   done
done
