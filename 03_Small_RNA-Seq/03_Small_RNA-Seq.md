1. Small RNA-Seq data:
Raw fastq files are accessible at NCBI's SRA under accession number SRP125901. Adapter-trimmed fastq files are available at:
https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/03_Small_RNA-Seq//Trimmed_fastq.tar
File nomenclature in that archive: STAGE_LIB_Trimmed.PF.R1.fastq.gz, with "STAGE" being "embryon_8h" for 8 hpf embryos (and idem for 15 hpf, 36 hpf and 60 hpf embryos), "femelle" being adult females, "male" being adult males; and "LIB" being either library #1, 2, 3 or 4 as described in Pinzón et al., 2018.

md5sum hash values for gzipped-compressed raw and trimmed fastq files are listed in 'md5sum_Seitz.txt'.

2. Genome mapping:
The B. lanceolatum genome assembly is available at:
https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/03_Small_RNA-Seq//Bl71nemr.fa.bz2
Bowtie2 index files for the Branchiostoma lanceolatum genome are available at:
https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/03_Small_RNA-Seq//B_lanceolatum_genome_bowtie2_index.tar.bz2

(download the adapter-trimmed fastq files and adapt the value of $PATH accordingly; download the bowtie2 index files and adapt the value of $PATH2 accordingly)

``for f in `ls $PATH/*Trimmed*`;do name=`echo $f | sed -e 's|.*/||' -e 's|_Trimmed.*||'`;bowtie2 -x $PATH2/B_lanceolatum -U $f --no-unal --quiet -S Mapping_$name'.sam';done``

3. Mapping on a database of B. lanceolatum abundant non-coding RNAs ('abundant_B_lanceolatum_ncRNAs.fa'; that database was assembled by blast-identification of orthologous loci for murine rRNAs, tRNAs, snRNAs, snoRNAs, scaRNAs and human 28S rRNA in the B. lanceolatum genome, then supplemented with sequences found after a sensitive homology search in the B. lanceolatum transcriptome, focusing on transcripts matched by high numbers of small RNA reads in our libraries):

``bowtie2-build abundant_B_lanceolatum_ncRNAs.fa abundant_B_lanceolatum_ncRNAs;for f in `ls $PATH/*Trimmed*`;do name=`echo $f | sed -e 's|.*/||' -e 's|_Trimmed.*||'`;bowtie2 -x abundant_B_lanceolatum_ncRNAs -U $f --no-unal --quiet -S ncRNA_mapping_$name'.sam';done``

4. Generating mapping statistics:

``echo "Sample Library Reads Trimmed_reads Genome-matching_reads ncRNA-matching_reads" > Statistics.dat;for sample in embryon_8h embryon_15h embryon_36h embryon_60h femelle male;do for lib in `seq 1 4`;do n1=`zcat $PATH/$sample'_'$lib'_S'*.fastq.gz | wc -l`;n1=`echo $n1"/4" | bc`;n2=`zcat $PATH/$sample'_'$lib'_Trimmed'*.fastq.gz | wc -l`;n2=`echo $n2"/4" | bc`;f=Mapping_$sample'_'$lib'.sam';n3=`sed '1,/^@PG\t/ d' $f | wc -l`;f=ncRNA_mapping_$sample'_'$lib'.sam';n4=`sed '1,/^@PG\t/ d' $f | wc -l`;echo $sample $lib $n1 $n2 $n3 $n4;done;done >> Statistics.dat``


5. Size distribution and logo analysis of genome-matching reads that do not match abundant non-coding RNAs:

``for sample in embryon_8h embryon_15h embryon_36h embryon_60h femelle male;do for lib in `seq 1 4`;do f=ncRNA_mapping_$sample'_'$lib'.sam';sed '1,/^@PG\t/ d' $f | awk '{print $1}' > ncRNA-matching_reads_$sample'_'$lib;f=Mapping_$sample'_'$lib'.sam';sed '1,/^@PG\t/ d' $f | awk '{print $1}' > genome-matching_reads_$sample'_'$lib;cat genome-matching_reads_$sample'_'$lib ncRNA-matching_reads_$sample'_'$lib ncRNA-matching_reads_$sample'_'$lib | sort | uniq -c | grep '^ *1 ' | awk '{print $2}' > genome-matching_reads_non_ncRNA-matching_reads_$sample'_'$lib;./Module_extract_from_fastq.pl `ls $PATH/$sample'_'$lib'_Trimmed'*` genome-matching_reads_non_ncRNA-matching_reads_$sample'_'$lib > genome-matching_reads_non_ncRNA-matching_reads_$sample'_'$lib'.fastq';done;done;for f in `ls genome-matching_reads_non_ncRNA-matching_reads_*.fastq`;do ./Module_converts_fastq_to_fasta.pl $f;done;MIN_LENGTH=18;MAX_LENGTH=30;for f in `ls genome-matching_reads_non_ncRNA-matching_reads_*.fa`;do name=`echo $f | sed -e 's|genome-matching_reads_non_ncRNA-matching_reads_||' -e 's|\.fa$||'`;echo "#Size Number_of_genome-matching_non_abundant_ncRNA-matching_reads" > 'Size_distribution_'$name'.dat';for size in `seq $MIN_LENGTH $MAX_LENGTH`;do grep -B 1 '^[ACGT]\{'$size'\}$' $f | grep -v '\-\-' > tmp_$size'_'$name'.fa';echo $size `grep -c '>' tmp_$size'_'$name'.fa'`;seqlogo -f tmp_$size'_'$name'.fa' -F EPS -o $size'-mer_logo_'$name -abcMnY;done >> 'Size_distribution_'$name'.dat';done``

6. Annotation of orthologous B. lanceolatum loci for known B. floridae or B. belcheri pre-miRNA hairpins, size distribution and logo analysis:

``./Fuses_lines_clean.pl hairpinMar18.fa | egrep -A 1 "^>bfl\-|^>bbe\-" | grep -v '\-\-' > Branchiostoma_known_hairpins.fa;blastn -db $PATH2/Bl71nemr.fa -query Branchiostoma_known_hairpins.fa -evalue 0.01 -word_size 4 -outfmt 6 > blast_output_hairpins.txt;./Module_extracts_miRNA_hairpins.pl blast_output_hairpins.txt``

Extracted sequences are in 'Orthologs_to_known_hairpins.fa'. Note that the dataset is redundant (if a lanceolatum sequence matches both a floridae and a belcheri sequence, it will appear twice in the list, maybe not with the exact same locus coordinates: there can be differences on the exact beginning and end).

``bowtie2-build Orthologs_to_known_hairpins.fa hairpins;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do for lib in `seq 1 4`;do bowtie2 -x hairpins -U genome-matching_reads_non_ncRNA-matching_reads_$sample'_'$lib'.fastq' --quiet --no-unal -S Hairpin-matching_reads_$sample'_'$lib'.sam';./Module_logo_for_hairpin-matching.pl Hairpin-matching_reads_$sample'_'$lib'.sam' > Size_dist_hairpin-matching_$sample'_'$lib'.dat';done;done``


7. Among the orthologous loci for known B. floridae or B. belcheri pre-miRNA hairpins, select those with classical pre-miRNA features (stable unbranched hairpins generating mostly 21-23-mers from their arms):

Unification of (possibly redundant) sequences in 'Orthologs_to_known_hairpins.fa':

``grep -v '^>' Orthologs_to_known_hairpins.fa | sort | uniq | nl | sed -e 's|^ *|>|' -e 's|\t|\
|' > tmp_unified.fa;makeblastdb -dbtype nucl -in tmp_unified.fa;blastn -db tmp_unified.fa -query tmp_unified.fa -evalue 0.01 -word_size 10 -outfmt 6 > blast_output_unification.txt;./Module_unification.pl blast_output_unification.txt tmp_unified.fa > unified1.fa;wc -l unified1.fa``

Result: 12582 lines in 'unified1.fa'.

``makeblastdb -dbtype nucl -in unified1.fa;blastn -db unified1.fa -query unified1.fa -evalue 0.01 -word_size 10 -outfmt 6 > blast_output_unification2.txt;./Module_unification.pl blast_output_unification2.txt unified1.fa > unified2.fa;wc -l unified2.fa``

Result: 12582 lines in 'unified2.fa' (unification has already converged; that was not obvious, because of possible transitive unification of long chains or related sequences).

``bowtie2-build unified1.fa unified_hairpins;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do lib=1;bowtie2 -x unified_hairpins -U genome-matching_reads_non_ncRNA-matching_reads_$sample'_'$lib'.fastq' --quiet --no-unal -S Unified_hairpin-matching_reads_$sample'_'$lib'.sam';done;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do lib=1;sed '1,/^@PG\t/ d' Unified_hairpin-matching_reads_$sample'_1.sam' | awk '$2==0 {print $3,$4,$10}' | sort | uniq -c > hairpin_mapping_$sample'_1.dat';done;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do echo "Unified_hairpin_ID bp_start read_seq count" > Hairpin_mapping_$sample'_1.dat';awk '{print $2,$3,$4,$1}' hairpin_mapping_$sample'_1.dat' >> Hairpin_mapping_$sample'_1.dat';rm -f hairpin_mapping_$sample'_1.dat';done;for id in `grep '^>' unified1.fa | sed -e 's|^> *||' -e 's| .*||'`;do seq=`grep -A 1 '^> *'$id' ' unified1.fa | tail -1`;l=`echo -n $seq | wc -m`;echo $id $l;done > Hairpin_lengths.dat;cat unified1.fa | RNAfold > Hairpin_folding.txt;sed -e '/^>/ s| .*||' -e '/^>/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' -e '/^>/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' -e 's|^>||' -e 's| ( *\(\-*[0-9\.]*\))$| \1|' Hairpin_folding.txt | tr acgu ACGU > Hairpin_structure.dat;R CMD BATCH R_commands_hairpin_read_profile``

These commands will have generated a series of 56 PDF files named 'Hairpin_*_read_coverage.pdf' (each one describes small RNA read coverage along the sequence of one of the 56 selected hairpins).

``for hairpin in `ls Hairpin_*.pdf | sed -e 's|^Hairpin_||' -e 's|_read_coverage.pdf$||'`;do grep -A 1 '^> *'$hairpin' ' unified1.fa;done > Convincing_pre-miRNA_hairpins.fa;bowtie2-build Convincing_pre-miRNA_hairpins.fa Convincing_hairpins;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do lib=1;bowtie2 -x Convincing_hairpins -U genome-matching_reads_non_ncRNA-matching_reads_$sample'_'$lib'.fastq' --quiet --no-unal -S Convincing_hairpin-matching_reads_$sample'_'$lib'.sam';done;./Script_miRNAs_in_development.sh``

These commands will have generated a series of 56 PDF files named 'Developmental_profile_miRNAs_from_hairpin_*.pdf' (each one shows the dynamics of 5´ arm and 3´ arm miRNA from one hairpin, during development). It will also have generated a TeX file named 'table_miRNAs_development.tex' for the preparation of Supplementary Table 1.


8. Size distribution and logo analysis of transcriptomic reads not matching abundant ncRNAs or pre-miRNA hairpins:

The B. lanceolatum transcriptome assembly is available at:
https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/03_Small_RNA-Seq//blanc_evm+_rn.fa.bz2

``for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do ./Script_transcriptome-matching_not_matching_ncRNA_or_hairpins.sh $sample;done;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do ./Chef3.sh $sample;done``


9. Size distribution and logo analysis of reads not matching abundant ncRNAs or pre-miRNA hairpins, while matching transcripts with long ORFs (Supplementary File S1 section 4):

Selection of transcripts with long ORFs (at least 100 codons, initiating on one of the first 3 AUG's):
``./Module_ORF_detector.pl Bla_annot_final_refTranscripts.fa``

(the resulting fasta file, with selected transcripts, is named 'Transcripts_with_convicing_ORFs_in_Bla_annot_final_refTranscripts.fa').

``bowtie2-build Transcripts_with_convicing_ORFs_in_Bla_annot_final_refTranscripts.fa Convincing_ORFs;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do ./Script_convincing_transcriptome-matching_not_matching_ncRNA_or_hairpins.sh $sample;done``


10. Size distribution and logo analysis of extragenomic, extratranscriptomic reads (Supplementary File S1 section 5):

``for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do for lib in `seq 1 4`;do bowtie2 -x /mnt/data/home/herve.seitz/Genomes/Bowtie2_indexes/Branchiostoma_lanceolatum/B_lanceolatum_transcriptome -U ../Extragenomic_reads_in_$sample'_'$lib'.fastq' --un Extragenomic_extratranscriptomic_reads_$sample'_'$lib'.fastq' > /dev/null;done;done
for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do for lib in `seq 1 4`;do ../Module_converts_fastq_to_fasta.pl Extragenomic_extratranscriptomic_reads_$sample'_'$lib'.fastq';grep -B 1 '^[ACGT]\{'18,30'\}$' Extragenomic_extratranscriptomic_reads_$sample'_'$lib'.fa' | grep -v '\-\-' > Extragenomic_extratranscriptomic_18_to_30nt_reads_$sample'_'$lib'.fa';done;done;for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do for lib in `seq 1 4`;do echo "Size Number_of_reads" > Size_distribution_extragenomic_extratranscriptomic_reads_$sample'_'$lib'.dat';for size in `seq 18 30`;do display=$size;grep -B 1 '^[ACGT]\{'$size'\}$' Extragenomic_extratranscriptomic_18_to_30nt_reads_$sample'_'$lib'.fa' | sed '/^>/ !s|T|U|g' | grep -v '\-\-' > $size'-mers_extragenomic_extratranscriptomic_reads_'$sample'_'$lib'.fa';nb_seq=`grep -c '>' $size'-mers_extragenomic_extratranscriptomic_reads_'$sample'_'$lib'.fa'`;display=`echo $display $nb_seq`;echo $display >> Size_distribution_extragenomic_extratranscriptomic_reads_$sample'_'$lib'.dat';seqlogo -f $size'-mers_extragenomic_extratranscriptomic_reads_'$sample'_'$lib'.fa' -F EPS -o Logo_$size'-mers_extragenomic_extratranscriptomic_reads_'$sample'_'$lib -abcMnY;done;done;done``


11. Analysis of C. elegans libraries (Supplementary File S1 section 6):

Download of sequence files for the three libraries (NCBI's GEO accession numbers GSM455391, GSM455392 and GSM455393):

``wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM455nnn/GSM455391/suppl/GSM455391%5FWT%5FTerminator%5Fseq%5FCCM%2Ds1%2Etar%2Egz;wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM455nnn/GSM455392/suppl/GSM455392%5FWT%5FTerminator%5Fseq%5FCCM%2Ds2%2Etar%2Egz;wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM455nnn/GSM455393/suppl/GSM455393%5FWT%5FTerminator%5Fseq%5FJCC%2Ds2%2Etar%2Egz``

Conversion into FASTA:

``for f in `ls GSM45539*.tar.gz1`;do tar -xzf $f;done;awk '{print $5}' s1.txt | grep -v '\.' | sort | uniq -c | nl | sed 's|^ *\([0-9]*\)\t *\([0-9]*\) \(.*\)|> seq_\1 (abundance=\2)\
\3|' > GSM455391.fa;awk '{print $5}' N2ind_2pm.txt | grep -v '\.' | sort | uniq -c | nl | sed 's|^ *\([0-9]*\)\t *\([0-9]*\) \(.*\)|> seq_\1 (abundance=\2)\
\3|' > GSM455392.fa;awk '{print $2}' s_2_eland_result.txt | grep -v '\.' | sort | uniq -c | nl | sed 's|^ *\([0-9]*\)\t *\([0-9]*\) \(.*\)|> seq_\1 (abundance=\2)\
\3|' > GSM455393.fa``

Adapter removal:

``for f in `ls GSM*.fa`;do ~/Annotation_pipeline/cutadapt-1.4.2/bin/cutadapt -a CTGTAGGCACCATCAAT $f -o trimmed_$f --discard-untrimmed;done;for f in `ls trimmed*.fa`;do grep -B 1 '^GGG' $f | grep -v '\-\-' | sed -e 's|^GGG||' -e 's|^> |>|' -e '/^>/ s| |_|g' > fully_$f;done # See Gu et al., 2009's online methods: "In  the  5´-ligation-independent  procedure,  Solexa  reads  containing  the  5´  linker  sequences (GGG) at position 1-3 and a perfect match to the first 6 nt of the 3´ linker (CTGTAG) were used to extract the inserted sequences" (so there's a 5´ GGG to be removed in addition to the 3´ adapter sequence)``

Assembly of a C. elegans "abundant non-coding RNA database" (on May 23, 2018):
a) Ribosomal RNAs: found by BLAST'ing U13369 on the C. elegans sequences, and with a search for "Caenorhabditis elegans [orgn] 5S" in NCBI's Nucleotide database.
b) Transfer RNAs: downloaded from the Genomic tRNA database (http://gtrnadb.ucsc.edu/genomes/eukaryota/Celeg_WS220/ce10-tRNAs.fa).
c) Small nuclear RNAs: found with a search for 'Caenorhabditis elegans [orgn] "small nuclear RNA"', in NCBI's Nucleotide database, setting sequence length between 1 and 1000 nt, then browsing for the names of spliceosomal snRNAs.
d) snoRNAs: downloaded from http://snoopy.med.miyazaki-u.ac.jp/snorna_db.cgi?mode=sno_search&organism=Caenorhabditis_elegans (selecting for "curated" sequences, then hand-checking each of their tick boxes).
e) Concatenation into a unique fasta file, named 'cel_abundant_ncRNAs.fa'.

C. elegans genome and transcriptome download (on June 1, 2017; "ce11" assembly), Bowtie2 index preparation:

``wget http://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/refMrna.fa.gz;gunzip refMrna.fa.gz;bowtie2-build refMrna.fa ce11_transcriptome;wget http://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz;tar -xzf chromFa.tar.gz;bowtie2-build chrI.fa,chrII.fa,chrIII.fa,chrIV.fa,chrM.fa,chrV.fa,chrX.fa ce11_genome``

pre-miRNA hairpin sequence download from miRBase (v.22, dated March 2018), Bowtie2 index preparation for C. elegans hairpins:

``wget ftp://mirbase.org/pub/mirbase/22/hairpin.fa.gz;gunzip hairpin.fa.gz;mv hairpin.fa hairpinMar18.fa;./Fuses_lines_clean.pl hairpinMar18.fa | grep -A 1 '^>cel\-' | grep -v '\-\-' | sed '/^>/ !s|U|T|g' > cel_hairpinMar18.fa;bowtie2-build cel_hairpinMar18.fa cel_pre-miRNA_hairpins``

Mapping the three Small RNA-Seq libraries on the C. elegans "abundant non-coding RNA database", on the genome and on pre-miRNA hairpins:

``bowtie2-build cel_abundant_ncRNAs.fa abundant_C_elegans_ncRNAs;for f in `ls fully_trimmed*.fa`;do bowtie2 -x abundant_C_elegans_ncRNAs -f -U $f -k 1 --no-unal -S abundant_cel_ncRNA_mapping_`echo $f | sed -e 's|fully_trimmed_||' -e 's|\.fa$||'`.sam;bowtie2 -x ce11_genome -f -U $f -k 1 --no-unal -S Genome_mapping_`echo $f | sed -e 's|fully_trimmed_||' -e 's|\.fa$||'`.sam;done;for f in `ls fully_trimmed*.fa`;do bowtie2 -x cel_pre-miRNA_hairpins -f -U $f -k 1 --no-unal -S hairpin_mapping_`echo $f | sed -e 's|fully_trimmed_||' -e 's|\.fa$||'`.sam;done``

Extraction of read IDs for genome-, abundant ncRNA- or pre-miRNA hairpin-matching reads:

``for f in `ls Genome_mapping_*.sam`;do sed '1,/^@PG\t/ d' $f | awk -F '\t' '{print $1}' > `echo $f | sed -e 's|^|Read_IDs_|' -e 's|\.sam$|.txt|'`;done;for f in `ls abundant_cel_ncRNA_mapping_*.sam`;do sed '1,/^@PG\t/ d' $f | awk -F '\t' '$2==0 {print $1}' > `echo $f | sed -e 's|^|Read_IDs_|' -e 's|\.sam$|.txt|'`;done;for f in `ls hairpin_mapping_*.sam`;do sed '1,/^@PG\t/ d' $f | awk -F '\t' '$2==0 {print $1}' > `echo $f | sed -e 's|^|Read_IDs_sense_|' -e 's|\.sam$|.txt|'`;sed '1,/^@PG\t/ d' $f | awk -F '\t' '$2==16 {print $1}' > `echo $f | sed -e 's|^|Read_IDs_antisense_|' -e 's|\.sam$|.txt|'`;done``

Preparation of fasta files for (genome-matching, not abundant ncRNA-matching), (pre-miRNA hairpin-matching), (not abundant ncRNA-matching, not pre-miRNA hairpin-matching) and (extragenomic) reads:

``for lib in GSM455391 GSM455392 GSM455393;do cat Read_IDs_Genome_mapping_$lib'.txt' Read_IDs_abundant_cel_ncRNA_mapping_$lib'.txt' Read_IDs_abundant_cel_ncRNA_mapping_$lib'.txt' | sort | uniq -c | grep '^ *1 ' | awk '{print $2}' > genomic_reads_not_matching_abundant_cel_ncRNA_$lib'.txt';grep '^>' fully_trimmed_$lib'.fa' | sed 's|^> *||' > all_reads_$lib'.txt';cat all_reads_$lib'.txt' Read_IDs_abundant_cel_ncRNA_mapping_$lib'.txt' Read_IDs_abundant_cel_ncRNA_mapping_$lib'.txt' Read_IDs_sense_hairpin_mapping_$lib'.txt' Read_IDs_sense_hairpin_mapping_$lib'.txt' | sort | uniq -c | grep '^ *1 ' | awk '{print $2}' > reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_$lib'.txt';cat all_reads_$lib'.txt' Read_IDs_Genome_mapping_$lib'.txt' Read_IDs_Genome_mapping_$lib'.txt' | sort | uniq -c | grep '^ *1 ' | awk '{print $2}' > extragenomic_reads_$lib'.txt';for f in genomic_reads_not_matching_abundant_cel_ncRNA_$lib'.txt' Read_IDs_sense_hairpin_mapping_$lib'.txt' Read_IDs_antisense_hairpin_mapping_$lib'.txt' reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_$lib'.txt' extragenomic_reads_$lib'.txt';do ./Module_extract_from_fasta.pl fully_trimmed_$lib'.fa' $f > `echo $f | sed 's|\.txt$|.fa|'`;done;done``

Transcriptome-matching reads that don't match abundant ncRNAs or pre-miRNA hairpins:

``for f in `ls reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_*'.fa'`;do bowtie2 -x ce11_transcriptome -f -U $f -k 1 --no-unal -S Transcriptome_mapping_`echo $f | sed -e 's|\.fa$||'`.sam;done``

Convincing ORF transcript-matching reads that don't match abundant ncRNAs or pre-miRNA hairpins:

``./Module_ORF_detector_lower_case.pl refMrna.fa;./Module_TXT_to_FA.pl Transcripts_with_convicing_ORFs_in_refMrna.txt refMrna.fa > Transcripts_with_convicing_ORFs_in_refMrna.fa;bowtie2-build Transcripts_with_convicing_ORFs_in_refMrna.fa cel_Convincing_ORFs;for f in `ls reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_*'.fa'`;do bowtie2 -x cel_Convincing_ORFs -f -U $f -k 1 --no-unal -S Convincing_ORF_transcript_mapping_`echo $f | sed -e 's|\.fa$||'`.sam;done``

Identification of (extragenomic and extratranscriptomic) and of (extragenomic and transcriptomic) reads:

``for f in `ls extragenomic_reads_*.fa`;do bowtie2 -x ce11_transcriptome -f -U $f -k 1 -S Transcriptome_mapping_`echo $f | sed -e 's|\.fa$||'`.sam;done``

Extraction of fasta files for (transcriptome-matching, not matching abundant ncRNAs or pre-miRNA hairpins), (convincing ORF transcript-matching, not matching abundant ncRNAs or pre-miRNA hairpins), (extragenomic and extratranscriptomic) and (extragenomic and transcriptomic) reads:

``for lib in GSM455391 GSM455392 GSM455393;do for analysis in Transcriptome_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins Convincing_ORF_transcript_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins;do sed '1,/^@PG\t/ d' $analysis'_'$lib'.sam' | awk '$2==0 {print ">"$1"\n"$10}' > sense_mapping_$analysis'_'$lib'.fa';sed '1,/^@PG\t/ d' $analysis'_'$lib'.sam' | awk '$2==16 {print ">"$1,$10}' > tmp_rev_$lib;awk '{print $1}' tmp_rev_$lib > tmp_column1_$lib;awk '{print $2}' tmp_rev_$lib | rev | tr ACGT TGCA > tmp_column2_$lib;paste -d '\n' tmp_column1_$lib tmp_column2_$lib > antisense_mapping_$analysis'_'$lib'.fa';done;sed '1,/^@PG\t/ d' Transcriptome_mapping_extragenomic_reads_$lib'.sam' | awk '$2==0 {print ">"$1"\n"$10}' > sense_mapping_transcriptomic_extragenomic_reads_$lib'.fa';sed '1,/^@PG\t/ d' Transcriptome_mapping_extragenomic_reads_$lib'.sam' | awk '$2==16 {print ">"$1,$10}' > tmp_rev_$lib;awk '{print $1}' tmp_rev_$lib > tmp_column1_$lib;awk '{print $2}' tmp_rev_$lib | rev | tr ACGT TGCA > tmp_column2_$lib;paste -d '\n' tmp_column1_$lib tmp_column2_$lib > antisense_mapping_transcriptomic_extragenomic_reads_$lib'.fa';sed '1,/^@PG\t/ d' Transcriptome_mapping_extragenomic_reads_$lib'.sam' | awk '$2==4 {print ">"$1"\n"$10}' > extratranscriptomic_extragenomic_reads_$lib'.fa';done``

Size distribution and logo analysis:

``for lib in GSM455391 GSM455392 GSM455393;do cp Read_IDs_sense_hairpin_mapping_$lib'.fa' sense_Read_IDs_hairpin_mapping_$lib'.fa';cp Read_IDs_antisense_hairpin_mapping_$lib'.fa' antisense_Read_IDs_hairpin_mapping_$lib'.fa';for file in genomic_reads_not_matching_abundant_cel_ncRNA_$lib'.fa' sense_Read_IDs_hairpin_mapping_$lib'.fa' antisense_Read_IDs_hairpin_mapping_$lib'.fa' sense_mapping_Transcriptome_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_$lib'.fa' antisense_mapping_Transcriptome_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_$lib'.fa' sense_mapping_Convincing_ORF_transcript_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_$lib'.fa' antisense_mapping_Convincing_ORF_transcript_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins_$lib'.fa' extratranscriptomic_extragenomic_reads_$lib'.fa' sense_mapping_transcriptomic_extragenomic_reads_$lib'.fa' antisense_mapping_transcriptomic_extragenomic_reads_$lib'.fa';do for size in `seq 18 30`;do grep -B 1 '^[ACGT]\{'$size'\}$' $file | grep -v '\-\-' | sed '{
N
s|>\(.*\)abundance=\([0-9]*\))\n\(.*\)|for i in \`seq 1 \2\`;do echo ">\1abundance=\2)";echo "\3";done|
}' > script.sh;chmod u+x script.sh;./script.sh | sed '/^>/ !s|T|U|g' > $size'-mers_in_'$file;seqlogo -f $size'-mers_in_'$file -F EPS -o Logo_$size'-mers_'`echo $file | sed 's|\.fa$||'` -abcMnY;done;done;for file in genomic_reads_not_matching_abundant_cel_ncRNA_$lib'.fa' extratranscriptomic_extragenomic_reads_$lib'.fa';do echo "Size Reads" > Size_distribution_`echo $file | sed 's|\.fa$|.dat|'`;for size in `seq 18 30`;do n=`grep -B 1 '^[ACGT]\{'$size'\}$' $file | grep '^>' | sed -e 's|.*abundance=||' -e 's|)$||' | awk '{s+=$1} END {print s}'`;if test "$n" = "";then n=0;fi;echo $size $n;done >> Size_distribution_`echo $file | sed 's|\.fa$|.dat|'`;done;for analysis in Read_IDs_hairpin_mapping mapping_Transcriptome_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins mapping_Convincing_ORF_transcript_mapping_reads_not_matching_abundant_cel_ncRNA_not_matching_pre-miRNA_hairpins mapping_transcriptomic_extragenomic_reads;do echo "Size Sense_matching Antisense_matching" > Size_distribution_$analysis'_'$lib'.dat';for size in `seq 18 30`;do s=`grep -B 1 '^[ACGT]\{'$size'\}$' sense_$analysis'_'$lib'.fa' | grep '^>' | sed -e 's|.*abundance=||' -e 's|)$||' | awk '{s+=$1} END {print s}'`;a=`grep -B 1 '^[ACGT]\{'$size'\}$' antisense_$analysis'_'$lib'.fa' | grep '^>' | sed -e 's|.*abundance=||' -e 's|)$||' | awk '{s+=$1} END {print s}'`;if test "$s" = "";then s=0;fi;if test "$a" = "";then a=0;fi;echo $size $s $a;done >> Size_distribution_$analysis'_'$lib'.dat';done;done;for f in `ls Logo*GSM*.eps`;do epstopdf $f;done;echo "Sample Genomic_reads_not_matching_abundant_cel_ncRNAs" > cel_Statistics.dat;for lib in GSM455391 GSM455392 GSM455393;do echo $lib `sed -e 's|.*abundance=||' -e 's|)$||' genomic_reads_not_matching_abundant_cel_ncRNA_$lib'.txt' | awk '{s+=$1} END {print s}'`;done >> cel_Statistics.dat;for lib in GSM455391 GSM455392 GSM455393;do for cmd in `ls R_commands_worm_*`;do Rscript $cmd $lib;done;done``


12. Analysis of extragenomic reads matching exon-exon junctions (Table 1 and Supplementary Figure 3):

Mapping introns and exons on the genome, using UCSC blat:

``blat ~/Genomes/B_lanceolatum/assembly/Bl71nemr.fa Bla_annot_final_refTranscripts.fa blat_output.psl``

(N.B.: the output file, 'blat_output.psl', can be downloaded in a bz2-compressed form from the address: https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/03_Small_RNA-Seq/blat_output.psl.bz2)

Then, select transcripts that contain at least one intron, and whose exon/intron structure is the same across all their genomic loci (some transcripts map at multiple loci, and their exon/intron boundaries are not consistent across loci: to be excluded), then select reads that map on exon-exon junctions:

``./Module_selects_from_blat_results.pl blat_output.psl;tail -n +2 Junction_and_exonic_reads_blat_output.dat | awk '{print $5}' | sort -g | uniq -c``

Results: 36 transcripts have reads mapping antisense to some of their exon-exon junctions. Transcripts with most such reads (when summing across all libraries and all developmental stages): BL05604_cuf0 (with 21 such reads; 20 of them are identical and the 21st is just a 1 nt variant of the same sequence; they all belong to "#3" libraries), BL00515_cuf0 (with 9 such reads; they are all identical, and they all belong to the "embryon_15h 1" library), BL16381_evm3 (with 6 such reads; 4 of them are identical, and the other 2 are 3´ extensions of the same sequence), and BL06097_cuf6 (with 4 such reads; 2 of them are identical, and the other 2 are 3´ extensions of the same sequence). In order to draw a map of these 4 loci (Supplementary Figure S3 only shows the top 2, for simplicity): generate data for genomic-matching reads:

``echo "Scaffold Strand Start End Exon_sizes Exon_starts Transcript" > RdRP_template_loci.txt;egrep -w "BL05604_cuf0|BL00515_cuf0|BL16381_evm3|BL06097_cuf6" blat_output.psl | awk '{print $14,$9,$16,$17,$19,$21,$10}' >> RdRP_template_loci.txt;tail -n +2 RdRP_template_loci.txt | awk '{print $1"\t"$3"\t"$4"\t"$7"\t1\t"$2}' > RdRP_template_loci.bed;sample=embryon_8h;nohup ./Script_mapping_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_$sample'.txt' &;sample=embryon_15h;nohup ./Script_mapping_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_$sample'.txt' &;sample=embryon_36h;nohup ./Script_mapping_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_$sample'.txt' &;sample=embryon_60h;nohup ./Script_mapping_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_$sample'.txt' &;sample=femelle;nohup ./Script_mapping_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_$sample'.txt' &;sample=male;nohup ./Script_mapping_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_$sample'.txt' &``

After completion of these 6 jobs: generate data for exon-exon junction-matching reads, to add to the map:

``for transcript in `awk '{print $4}' RdRP_template_loci.bed`;do l=`grep -A 1 '^>'$transcript'$' Fused_Bla_annot_final_refTranscripts.fa | tail -1 | wc -m`;l=`echo $l"-1" | bc`;echo $transcript 1 $l;done | sed 's| |\t|g' > RdRP_template_transcripts.bed;sample=embryon_8h;nohup ./Script_mapping_extragenomic_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_extragenomic_$sample'.txt' &;sample=embryon_15h;nohup ./Script_mapping_extragenomic_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_extragenomic_$sample'.txt' &;sample=embryon_36h;nohup ./Script_mapping_extragenomic_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_extragenomic_$sample'.txt' &;sample=embryon_60h;nohup ./Script_mapping_extragenomic_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_extragenomic_$sample'.txt' &;sample=femelle;nohup ./Script_mapping_extragenomic_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_extragenomic_$sample'.txt' &;sample=male;nohup ./Script_mapping_extragenomic_on_RdRP_template_loci.sh $sample > nohup_RdRP_template_loci_extragenomic_$sample'.txt' &``

To assess how common sense and antisense exon-exon junction reads are (for Table 1):

``total_depth=`tail -n +2 ../Statistics.dat | awk '{s+=$5-$6} END {print s}'``
(this is the total number of genome-matching reads that do not match on abundant ncRNAs (for each of the 24 libraries); note that it includes read shorter than 18 nt or longer than 30 nt, like I have done for all my measurements of sequencing depth)
``echo "mRNA Number_of_sense_junction_reads(ppm) Number_of_antisense_junction_reads(ppm)" > Top_sense_junction_mapped_mRNAs.dat;tail -n +2 Junction_and_exonic_reads_blat_output.dat | awk '{print $3*1000000/'$total_depth',$5*1000000/'$total_depth',$1}' | sort -gr | awk '{print $3,$1,$2}' | head >> Top_sense_junction_mapped_mRNAs.dat;echo "mRNA Number_of_sense_junction_reads(ppm) Number_of_antisense_junction_reads(ppm)" > Top_antisense_junction_mapped_mRNAs.dat;tail -n +2 Junction_and_exonic_reads_blat_output.dat | awk '$5>0 {print $5*1000000/'$total_depth',$3*1000000/'$total_depth',$1}' | sort -gr | awk '{print $3,$2,$1}' | head >> Top_antisense_junction_mapped_mRNAs.dat``

Statistics transcriptome-wide (raw reads and ppm):

``awk '{s+=$3} END {print s}' Junction_and_exonic_reads_blat_output.dat # Number of raw reads mapping sense to junctions (result: 5683)``
``awk '{s+=$5} END {print s}' Junction_and_exonic_reads_blat_output.dat # Number of raw reads mapping antisense to junctions (result: 80)``
``awk '{s+=$3*1000000/'$total_depth'} END {print s}' Junction_and_exonic_reads_blat_output.dat`` (number of depth-normalized reads mapping sense to junctions (result: 94.2325 ppm))
``awk '{s+=$5*1000000/'$total_depth'} END {print s}' Junction_and_exonic_reads_blat_output.dat`` (number of depth-normalized reads mapping antisense to junctions (result: 1.32652 ppm))
