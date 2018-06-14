1. Identifying transcripts with most antisense exon-exon junction reads (relatively to sense exon-exon junction reads):

Using files 'Transcriptome_mapping_of_extragenomic_reads_*' from the analyses described in '03_Small_RNA-Seq.md'; these files are copied here (in the archive 'Transcriptome_mapping_of_extragenomic_reads.tar.bz2') for convenience.

``for sample in embryon_15h embryon_36h embryon_60h embryon_8h femelle male;do for lib in `seq 1 4`;do sed '1,/^@PG\t/ d' Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam' | awk '$2==0 {print $3,$10}' | grep ' [ACGTN]\{'18,30'\}$' > sense_junction_reads_$sample'_'$lib'.txt';sed '1,/^@PG\t/ d' Transcriptome_mapping_of_extragenomic_reads_$sample'_'$lib'.sam' | awk '$2==16 {print $3,$10}' | grep ' [ACGTN]\{'18,30'\}$' > antisense_junction_reads_$sample'_'$lib'.txt';echo "Transcript Sense_junction_reads Antisense_junction_reads" > Sense_vs_antisense_junction_reads_$sample'_'$lib'.dat';for g in `cat sense_junction_reads_$sample'_'$lib'.txt' antisense_junction_reads_$sample'_'$lib'.txt' | awk '{print $1}' | sort | uniq`;do s=`grep -c '^'$g' ' sense_junction_reads_$sample'_'$lib'.txt'`;a=`grep -c '^'$g' ' antisense_junction_reads_$sample'_'$lib'.txt'`;if test "$s" = "";then s=NA;fi;if test "$a" = "";then a=NA;fi;echo $g $s $a >> Sense_vs_antisense_junction_reads_$sample'_'$lib'.dat';done;R CMD BATCH R_ratio_antisense_sense_across_transcriptsi;sed -e 's|^\[1\] ||' -e 's|"||g' -e 's| |\
|g' Transcripts_enriched_in_antisense_junction_reads.txt > Prefered_RdRP_templates.txt``

The resulting file ('Prefered_RdRP_templates.txt') contains a list of 4 transcripts with particularly high antisense/sense fragment ratios. These are the 4 identified 'RdRP templates' shown in Figure 6.


2. mRNA abundance measurement using vast-tools (installed from https://github.com/vastgroup/vast-tools) for the adult RNA-Seq data (NCBI BioSample accession #SAMN09381006 and SAMN09381007), using the B. lanceolatum transcriptome database (which is available at: https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/03_Small_RNA-Seq//B_lanceolatum_transcriptome_bowtie2_index.tar.bz2):

``for sample in Ax1 Ax2 Ax7 Ax8;do vast-tools align `ls $sample'_'*_R1.fastq.gz` `ls $sample'_'*_R2.fastq.gz` -sp Bla --exprONLY;done``

The resulting files are copied here for convenience (they are named 'Ax1_ATCACG_L003_R1.cRPKM', 'Ax2_CGATGT_L003_R1.cRPKM', 'Ax7_CAGATC_L003_R1.cRPKM' and 'Ax8_ACTTGA_L003_R1.cRPKM'). For the joint analysis of these adult RNA-Seq datasets with embryonic and juvenile RNA-Seq data (from Marlétaz et al., submitted; these data are copied here for convenience, in file 'Bla_cRPKMs_only-49.tab.bz2'): verification that the genes are listed in the same order in each file (that will speed up downstream analyses):

``tail -n +2 Bla_cRPKMs_only-49.tab | awk '{print $1}' | md5sum;for f in `ls vast_out/expr_out/Ax*.cRPKM`;do awk '{print $1}' $f | md5sum;done``

(yes, they are)

Plotting gene expression dynamics for RdRP genes and the 4 identified RdRP templates:

``R CMD BATCH R_commands_R_commands_ANOVA_for_time_dependence_with_temporal_significativity``

That script will have generated PDF files showing the expression dynamics of the 6 RdRP genes and the 4 RdRP templates during development (using an arbitrary age for premetamorphic larvae, adult males and females: to be adapted by hand to make it explicit in the final figure).
