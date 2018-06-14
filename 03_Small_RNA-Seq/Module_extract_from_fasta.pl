#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter input file names (e.g., ./Module_extract_from_fasta.pl fully_trimmed_GSM455391.fa genomic_reads_not_matching_abundant_cel_ncRNA_GSM455391.txt).\n";
}
else
{
    $n=1;
    open(FASTA,$ARGV[0]);
    while(<FASTA>)
    {
	++$n;
	chomp;
	if ($n == 2)
	{
	    s/^> *//;
	    s/ .*//;
	    $name=$_;
	    $n=0;
	}
        if ($n == 1)
        {
            $l1{$name}=$_;
        }
    }
    close(FASTA);

    open(LIST,$ARGV[1]);
    while (<LIST>)
    {
	chomp;
	print ">$_\n$l1{$_}\n";
    }
    close(LIST);
}
