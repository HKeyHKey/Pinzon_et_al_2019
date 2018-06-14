#!/usr/bin/perl

if ($ARGV[0] eq '')
{
    print "Please enter input file name (e.g., ./Module_extracts_fastq_read_identifiers.pl trimmed_bam86.fastq).\n";
}
else
{
    $n=3;
    open(FASTQ,$ARGV[0]);
    while(<FASTQ>)
    {
	++$n;
	if ($n == 4)
	{
	    chomp;
	    s/^\@//;
	    s/ .*//;
	    print "$_\n";
	    $n=0;
	}
    }
    close(FASTQ);
}
