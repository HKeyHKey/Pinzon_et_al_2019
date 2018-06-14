echo "\\tableofcontents" > body.tex
echo "\\newpage" >> body.tex
for analysis in non_ncRNA_genomic hairpin-matching transcriptomic_not_ncRNA_not_hairpin-matching convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching extragenomic_extratranscriptomic transcriptomic_extragenomic
do case "$analysis" in "non_ncRNA_genomic") title="Genomic reads not matching abundant ncRNAs";;
		       "hairpin-matching") title="pre-miRNA hairpin-matching reads";;
                       "transcriptomic_not_ncRNA_not_hairpin-matching") title="Transcriptome-matching reads (excuding pre-miRNA and abundant ncRNA-matching reads)";;
                       "convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching") title="Reads matching RNAs with long ORF's";;
                       "extragenomic_extratranscriptomic") title="Extragenomic and extratranscriptomic reads";;
		       "transcriptomic_extragenomic") title="Extragenomic reads matching the transcriptome";;
   esac
   echo "\\section{"$title"}"
   for lib in `seq 1 4`
   do case "$lib" in "1") lib_title="total 5\'{} monophosphorylated small RNAs";;
                     "2") lib_title="3\'{} modified, 5\'{} monophosphorylated small RNAs";;
                     "3") lib_title="total 5\'{} hydroxyl or polyphosphorylated small RNAs";;
                     "4") lib_title="3\'{} modified, 5\'{} hydroxyl or polyphosphorylated small RNAs";;
      esac
      echo "\\subsection{Libraries \\#"$lib" ("$lib_title")}"
      for sample in embryon_8h embryon_15h embryon_36h embryon_60h femelle male
      do echo `echo $sample | sed -e 's|_| |g' -e 's|embryon|Embryo|' -e 's|male|Adult male|' -e 's|femelle|Adult female|'`", library "$lib":"
         echo ""
         echo ""
         echo "\\begin{tabular}{ccc}"
         size=18
         if test "$analysis" = "extragenomic_extratranscriptomic" -o "$analysis" = "non_ncRNA_genomic"
         then echo "&&18-mers:\\\\"
              echo "\\includegraphics[width=6cm]{Size_distribution_"$analysis"_reads_"$sample"_"$lib".pdf}&&\\includegraphics[width=6cm]{Logo_"$size"-mers_"$analysis"_reads_"$sample"_"$lib".pdf}\\\\"
              for i in `seq 1 4`
              do size1=`echo "16+3*"$i | bc`
                 size2=`echo "17+3*"$i | bc`
                 size3=`echo "18+3*"$i | bc`
                 echo $size1"-mers:&"$size2"-mers:&"$size3"-mers:\\\\"
                 echo "\\includegraphics[width=6cm]{Logo_"$size1"-mers_"$analysis"_reads_"$sample"_"$lib".pdf}&\\includegraphics[width=6cm]{Logo_"$size2"-mers_"$analysis"_reads_"$sample"_"$lib".pdf}&\\includegraphics[width=6cm]{Logo_"$size3"-mers_"$analysis"_reads_"$sample"_"$lib".pdf}\\\\"
              done
              echo "\\end{tabular}"
         else echo "&Sense reads:&18-mers:\\\\"
              echo "\\includegraphics[width=6cm]{Size_distribution_"$analysis"_reads_"$sample"_"$lib".pdf}&&\\includegraphics[width=6cm]{Logo_"$size"-mers_sense_"$analysis"_reads_"$sample"_"$lib".pdf}\\\\"
              for i in `seq 1 4`
              do size1=`echo "16+3*"$i | bc`
                 size2=`echo "17+3*"$i | bc`
                 size3=`echo "18+3*"$i | bc`
                 echo $size1"-mers:&"$size2"-mers:&"$size3"-mers:\\\\"
                 echo "\\includegraphics[width=6cm]{Logo_"$size1"-mers_sense_"$analysis"_reads_"$sample"_"$lib".pdf}&\\includegraphics[width=6cm]{Logo_"$size2"-mers_sense_"$analysis"_reads_"$sample"_"$lib".pdf}&\\includegraphics[width=6cm]{Logo_"$size3"-mers_sense_"$analysis"_reads_"$sample"_"$lib".pdf}\\\\"
              done
              echo "\\end{tabular}"
              echo "\\newpage"
              echo "\\begin{tabular}{ccc}"
              size=18
              echo "&Antisense reads:&18-mers:\\\\"
              echo "&&\\includegraphics[width=6cm]{Logo_"$size"-mers_sense_"$analysis"_reads_"$sample"_"$lib".pdf}\\\\"
              for i in `seq 1 4`
              do size1=`echo "16+3*"$i | bc`
                 size2=`echo "17+3*"$i | bc`
                 size3=`echo "18+3*"$i | bc`
                 echo $size1"-mers:&"$size2"-mers:&"$size3"-mers:\\\\"
                 echo "\\includegraphics[width=6cm]{Logo_"$size1"-mers_antisense_"$analysis"_reads_"$sample"_"$lib".pdf}&\\includegraphics[width=6cm]{Logo_"$size2"-mers_antisense_"$analysis"_reads_"$sample"_"$lib".pdf}&\\includegraphics[width=6cm]{Logo_"$size3"-mers_antisense_"$analysis"_reads_"$sample"_"$lib".pdf}\\\\"
              done
              echo "\\end{tabular}"
         fi
         echo "\\newpage"
      done
   done
done >> body.tex
# Below: needs to make sure no non-existing logo plot is invocated:
for size in `seq 18 30`;do for lib in `seq 1 4`;do for sample in embryon_8h embryon_15h embryon_36h embryon_60h femelle male;do if [ ! -f Logo_$size'-mers_sense_hairpin-matching_reads_'$sample'_'$lib'.pdf' ];then echo $size'-mers_sense_hairpin-matching_reads_'$sample'_'$lib;fi;if [ ! -f Logo_$size'-mers_antisense_hairpin-matching_reads_'$sample'_'$lib'.pdf' ];then echo $size'-mers_antisense_hairpin-matching_reads_'$sample'_'$lib;fi;done;done;done > tmp_missing
for size in `seq 18 30`;do for lib in `seq 1 4`;do for sample in embryon_8h embryon_15h embryon_36h embryon_60h femelle male;do if [ ! -f Logo_$size'-mers_sense_convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads_'$sample'_'$lib'.pdf' ];then echo $size'-mers_sense_convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads_'$sample'_'$lib;fi;if [ ! -f Logo_$size'-mers_antisense_convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads_'$sample'_'$lib'.pdf' ];then echo $size'-mers_antisense_convincing_ORF_transcriptomic_not_ncRNA_not_hairpin-matching_reads_'$sample'_'$lib;fi;done;done;done >> tmp_missing
for size in `seq 18 30`;do for lib in `seq 1 4`;do for sample in embryon_8h embryon_15h embryon_36h embryon_60h femelle male;do if [ ! -f Logo_$size'-mers_sense_transcriptomic_extragenomic_reads_'$sample'_'$lib'.pdf' ];then echo $size'-mers_sense_transcriptomic_extragenomic_reads_'$sample'_'$lib;fi;if [ ! -f Logo_$size'-mers_antisense_transcriptomic_extragenomic_reads_'$sample'_'$lib'.pdf' ];then echo $size'-mers_antisense_transcriptomic_extragenomic_reads_'$sample'_'$lib;fi;done;done;done >> tmp_missing
for f in `cat tmp_missing`;do sed -i 's|\\includegraphics\[width=6cm\]{Logo_'$f'\.pdf}|(no read)|g' body.tex;done
# Above: to make sure no non-existing logo plot is invocated
cat header.tex body.tex tailer.tex > Size_distributions_and_logos.tex
