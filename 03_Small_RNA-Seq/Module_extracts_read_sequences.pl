#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter input file names (a list of IDs and a fastq file containing these IDs; e.g., ./Module_extracts_read_sequences.pl dme_genome_matching_reads.txt trimmed_bam86.fastq).\n";
}
else
{
    $n=3;
    open(FASTQ,$ARGV[1]);
    while(<FASTQ>)
    {
	++$n;
	if ($n == 4)
	{
	    chomp;
	    s/^\@//;
	    s/ .*//;
	    $ID=$_;
	    $n=0;
	}
	if ($n == 1)
	{
	    chomp;
	    $seq{$ID}=$_;
	}
    }
    close(FASTQ);
   
#    $out_file=$ARGV[0];
#    $out_file=~s/\.txt$/.fa/;
#    open(OUT,">$out_file");
    open(LIST,$ARGV[0]);
    while(<LIST>)
    {
	chomp;
	if ($seq{$_} ne '')
	{
#	    print OUT ">$_\n$seq{$_}\n";
	    print ">$_\n$seq{$_}\n";
	}
    }
    close(LIST);
#    close(OUT);
}
