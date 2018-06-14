1. HMMer profile download from PFAM (on December 8, 2017):

``wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz;wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz``


2. Extraction of HMMer profiles for every annotated RdRP class:

``for acc in PF04197 PF00998 PF00972 PF02123 PF05919 PF00680 PF00978 PF17501 PF12426 PF08467 PF05788 PF03431 PF05183 PF00946 PF04196 PF07925 PF00603 PF00602 PF00604 PF03035;do sed '/^HMMER3\/f / {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' Pfam-A.hmm | sed '/^HMMER3\/f / {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed -e '/^HMMER3\/f .* ACC *'$acc'\./,/^\/\// !d' -e '/^HMMER3\/f / s| ACC|\
ACC|' -e '/^HMMER3\/f / s| NAME|\
NAME|' > $acc'.hmm'
done``


3. HMMer search with RdRP profiles on the 538 animal proteomes:

``for species in `ls Animal_proteomes/*.fa`;do name=`echo $species | sed -e 's|.*/||' -e 's|sta$||' -e 's|\.fa$||'`;for profile in `ls PF*.hmm`;do prof=`echo $profile | sed 's|\.hmm$||'`;hmmsearch $profile $species > $prof'_on_'$name'.out';done;done;for f in `ls PF*.out`;do bait=`echo $f | sed 's|_on_.*||'`;species=`echo $f | sed -e 's|.*_on_||' -e 's|\.out$||'`;for prot in `sed '/^Scores for complete sequences (score includes all domains):$/,/^ *$/ !d' $f | sed '1,/^ *-------  *------  *-----  *-------  *------  *-----  *----  *--  *--------  *-----------$/ d' | grep -v '\-\-\-\-\-\- inclusion threshold \-\-\-\-\-\-' | grep -v '^ *$' | awk '{print $9}'`;do grep -w -A 1 '^>'$prot Search_for_other_proteins/Fused_$species'.fa';done | sed -e 's|\t| |g' -e '/^>/ s|$| '$species'|' > Sequences_candidates_$bait'_'$species'.fa';done;cat Sequences_candidates_PF*.fa > Sequences_candidates_every_PFAM_RdRP.fa``

4. Screening for RdRP candidates with complete RdRP domains, using NCBI's "Conserved Domain" database:
Submission of 'Sequences_candidates_every_PFAM_RdRP.fa' to http//www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi on December 8, 2017. Results saved as 'Every_PFAM_RdRP_evening_hitdata_concise.txt', 'Every_PFAM_RdRP_evening_hitdata_standard.txt' and 'Every_PFAM_RdRP_evening_hitdata_full.txt'.

``for p in `ls PF*.hmm | sed -e 's|\.hmm$||' -e 's|PF|pfam|'`;do echo "*** "$p" ***";awk -F '\t' '$8=="'$p'" && $10==" - " {print $1}' Every_PFAM_RdRP_evening_hitdata_full.txt | sed -e 's|.* ||' -e 's|_proteome$||' | sort | uniq;echo "";done``

5. Displaying the taxonomy of RdRP-containing proteomes:

``for p in `ls PF*.hmm | sed -e 's|\.hmm$||' -e 's|PF|pfam|'`;do awk -F '\t' '$8=="'$p'" && $10==" - " {print $1}' Every_PFAM_RdRP_evening_hitdata_full.txt | sort | uniq > For_$p'_candidates_with_complete_domain.dat';done;for p in `ls PF*.hmm | sed -e 's|\.hmm$||' -e 's|PF|pfam|'`;do if [ -s For_$p'_candidates_with_complete_domain.dat' ];then nb=`cat For_$p'_candidates_with_complete_domain.dat' | wc -l`;for i in `seq 1 $nb`;do prot=`head -$i For_$p'_candidates_with_complete_domain.dat' | tail -1 | sed -e 's|^Q#[0-9]* - ||' -e 's| .*||'`;species=`head -$i For_$p'_candidates_with_complete_domain.dat' | tail -1 | sed -e 's|.* ||' -e 's|_proteome$||'`;if test "$species" = "floridae]";then species=Branchiostoma_floridae;fi;if [ ! -f Search_for_other_proteins/Fused_$species'_proteome.fa' ];then echo "Problem! I cannot find Fused_"$species"_proteome.fa";fi;grep -w -A 1 '^'$prot Search_for_other_proteins/Fused_$species'_proteome.fa' >> For_$p'_whole_sequence_candidates_'$species'.fa';done;fi;done;for clade in Brachiopoda Rhabditida Diplogasterida Spirurida Tylenchida Ascaridida Enoplea Oxyurida Scalidophora Tardigrada Merostomata Arachnida Myriapoda Crustacea Collembola Diptera Hymenoptera Amphiesmenoptera Coleoptera Hemiptera Polyneoptera Phthiraptera Gastropoda Bivalvia Cephalopoda Clitellata Polychaeta Sauropsida Mammalia Amphibia Coelacanthiformes Actinopterygii Chondrichthyes Hyperoartia Cephalochordata Tunicata Hemichordata Echinodermata Hydrozoa Anthozoa Myxozoa Porifera Mesozoa Placozoa Cestoda Trematoda Monogenea Rhabditophora;do echo "*** "$clade" ***";for species in `ls Animal_proteomes/*_proteome.fa | sed -e 's|Animal_proteomes/||' -e 's|_proteome.fa||'`;do if test `grep '^'$species': ' Animal_proteomes/Taxonomy_info | grep -wc $clade` -ne 0;then string=$species;for p in `ls PF*.hmm | sed -e 's|\.hmm$||' -e 's|PF|pfam|'`;do if [ -f For_$p'_whole_sequence_candidates_'$species'.fa' ];then n=`grep -c '^>' For_$p'_whole_sequence_candidates_'$species'.fa'`;string=`echo $string $p':'$n`;fi;done;echo $string;fi;done;done > taxo_out``

To plot the results described in 'taxo_out' as piecharts for Figure 1:

``for clade in `grep '\*' taxo_out | awk '{print $2}'`;do sed '/\*\*\* '$clade' \*\*\*/,/\*/ !d' taxo_out | grep -v '\*' > tmp;nb_species=`cat tmp | wc -l`;string=`echo $clade $nb_species`;sed -e 's|^[A-Z][a-z_]*||' -e 's|^ *||' -e 's|:[0-9]*||g' tmp | grep -v '^$' | sed 's| |\
|g' | sort | uniq > list_types;for t in `cat list_types`;do add=`grep -c ' '$t':' tmp`;string=`echo $string $t":"$add`;done;echo $string;done > Distribution_of_RdRP_types.dat;for clade in `awk '{print $1}' Distribution_of_RdRP_types.dat`;do line=`grep '^'$clade' ' Distribution_of_RdRP_types.dat`;total=`echo $line | awk '{print $2}'`;list=`echo $line | sed -e 's|^'$clade' '$total' *||' -e 's|:[0-9]*||g'`;for t in `echo $list`;do n=`echo $line | sed -e 's|.* '$t':||' -e 's| .*||'`;echo $clade $t $n"/"$total;done;done > Detailed_distribution_of_RdRP_types.dat;R CMD BATCH R_commands_piecharts``

6. To control for the effect of proteome quality (Suppl. Figure 1): exclusion of proteomes with suspiciously low numbers of long proteins (either less than 5,000 proteins of at least 500 amino acids; or less than 1,000 proteins of at least 1,000 amino acids):

``for clade in Brachiopoda Rhabditida Diplogasterida Spirurida Tylenchida Ascaridida Enoplea Oxyurida Scalidophora Tardigrada Merostomata Arachnida Myriapoda Crustacea Collembola Diptera Hymenoptera Amphiesmenoptera Coleoptera Hemiptera Polyneoptera Phthiraptera Gastropoda Bivalvia Cephalopoda Clitellata Polychaeta Sauropsida Mammalia Amphibia Coelacanthiformes Actinopterygii Chondrichthyes Hyperoartia Cephalochordata Tunicata Hemichordata Echinodermata Hydrozoa Anthozoa Myxozoa Porifera Mesozoa Placozoa Cestoda Trematoda Monogenea Rhabditophora;do echo "*** "$clade" ***";for species in `ls Animal_proteomes/*_proteome.fa | sed -e 's|Animal_proteomes/||' -e 's|_proteome.fa||'`;do if test `grep '^'$species': ' Animal_proteomes/Taxonomy_info | grep -wc $clade` -ne 0 -a `grep -wc $species Animal_proteomes/Species_with_at_least_1000_proteins_of_more_than_1000aa.txt` -ne 0;then string=$species;for p in `ls PF*.hmm | sed -e 's|\.hmm$||' -e 's|PF|pfam|'`;do if [ -f For_$p'_whole_sequence_candidates_'$species'.fa' ];then n=`grep -c '^>' For_$p'_whole_sequence_candidates_'$species'.fa'`;string=`echo $string $p':'$n`;fi;done;echo $string;fi;done;done > taxo_1000_of_1000_out``

``for clade in Brachiopoda Rhabditida Diplogasterida Spirurida Tylenchida Ascaridida Enoplea Oxyurida Scalidophora Tardigrada Merostomata Arachnida Myriapoda Crustacea Collembola Diptera Hymenoptera Amphiesmenoptera Coleoptera Hemiptera Polyneoptera Phthiraptera Gastropoda Bivalvia Cephalopoda Clitellata Polychaeta Sauropsida Mammalia Amphibia Coelacanthiformes Actinopterygii Chondrichthyes Hyperoartia Cephalochordata Tunicata Hemichordata Echinodermata Hydrozoa Anthozoa Myxozoa Porifera Mesozoa Placozoa Cestoda Trematoda Monogenea Rhabditophora;do echo "*** "$clade" ***";for species in `ls Animal_proteomes/*_proteome.fa | sed -e 's|Animal_proteomes/||' -e 's|_proteome.fa||'`;do if test `grep '^'$species': ' Animal_proteomes/Taxonomy_info | grep -wc $clade` -ne 0 -a `grep -wc $species Animal_proteomes/Species_with_at_least_5000_proteins_of_more_than_500aa.txt` -ne 0;then string=$species;for p in `ls PF*.hmm | sed -e 's|\.hmm$||' -e 's|PF|pfam|'`;do if [ -f For_$p'_whole_sequence_candidates_'$species'.fa' ];then n=`grep -c '^>' For_$p'_whole_sequence_candidates_'$species'.fa'`;string=`echo $string $p':'$n`;fi;done;echo $string;fi;done;done > taxo_5000_of_500_out``


Plotting piecharts for Supplementary Figure 1 (similar to those of Figure 1A, but after exclusion of proteomes with low numbers of long proteins):

``for select in 1000_of_1000 5000_of_500;do for clade in `grep '\*' taxo_$select'_out' | awk '{print $2}'`;do sed '/\*\*\* '$clade' \*\*\*/,/\*/ !d' taxo_$select'_out' | grep -v '\*' > tmp;nb_species=`cat tmp | wc -l`;string=`echo $clade $nb_species`;sed -e 's|^[A-Z][a-z_]*||' -e 's|^ *||' -e 's|:[0-9]*||g' tmp | grep -v '^$' | sed 's| |\
|g' | sort | uniq > list_types;for t in `cat list_types`;do add=`grep -c ' '$t':' tmp`;string=`echo $string $t":"$add`;done;echo $string;done > Distribution_of_RdRP_types_$select'.dat';for clade in `awk '{print $1}' Distribution_of_RdRP_types_$select'.dat'`;do line=`grep '^'$clade' ' Distribution_of_RdRP_types_$select'.dat'`;total=`echo $line | awk '{print $2}'`;list=`echo $line | sed -e 's|^'$clade' '$total' *||' -e 's|:[0-9]*||g'`;for t in `echo $list`;do n=`echo $line | sed -e 's|.* '$t':||' -e 's| .*||'`;echo $clade $t $n"/"$total;done;done > Detailed_distribution_of_RdRP_types_$select'.dat';done``

then plotting with R, using the commands described in 'R_commands_piecharts'.

7. Re-sequencing of the Branchiostoma lanceolatum BL09945 locus:

Amplification of the BL09945 locus with oligos GCCTTATTCGTCTATGGCTGTT and CAGCTACGCCTGTGTTTACG, cloning in pGEM-T easy then sequencing several clones of the resulting plasmids with the following oligos:

a) with T7 promoter-specific (TAATACGACTCACTATAGGG) and SP6 promoter-specific (CATTTAGGTGACACTATAG) oligos: results are archived in 'Order_11104254864-1_Results.zip';

b) with a series of internal oligos: d1835 (GGAAGGCACGCAGTTCTTCGG), d1836 (CTGGCTGTTTAGCTGTGGAAG), d1837 (CCACTTAGGCTGCCTTTGGGC), d1838 (ACTGACTGTTCAACTGTCGGA), d1839 (CAAGATGCGCGTGTGCAATCT), d1840 (GCCTTTGGGCAAACTGTCGCG) and d1841 (CAAGGAAGAGCGATCAAAGGT): results are archived in 'Order_11104283914-1_Results.zip'.

Sequences of clones #4 and #5 were assembled using pregap4 and gap4 (with default parameters) and the resulting fasta files are named 'Clone4.fa' and 'Clone5.fa' respectively. Removal of flanking vector sequences:

``~/Fuses_lines_clean.pl Sequencing_23012018/Clone4.fa | grep -A 1 SP6 | tr acgt ACGT | sed -e 's|.*CCCATATGGTCGACCTGCAGGCGGCCGCGAATTCACTAGTGATTC*A*||' -e 's|ATCGAATTCCCGCGGCCGCCATGG.*||' -e '/^>/ s|.*|>Variant_1 (clone 4)|' > For_GenBank_submission.fa;~/Fuses_lines_clean.pl Sequencing_23012018/Clone5.fa | grep -A 1 SP6 | tr acgt ACGT | sed -e 's|.*CCCATATGGTCGACCTGCAGGCGGCCGCGAATTCACTAGTGATTC*A*||' -e 's|ATCGAATTCCCGCGGCCGCCATGG.*||' -e '/^>/ s|.*|>Variant_2 (clone 5)|' >> For_GenBank_submission.fa``

An artifactual 1-nt deletion in both clones around bp 1180 (ATCGAAGATTTCATT was changed into ATCGAAGATTTCAT) was corrected by hand after chromatogram inspection:

``sed -i 's|ATCGAAGATTTCAT|ATCGAAGATTTCATT|' For_GenBank_submission.fa``
