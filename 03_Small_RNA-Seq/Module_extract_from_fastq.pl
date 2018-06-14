#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter input file names (e.g., ./Module_extract_from_fastq.pl trimmed_TM3hs20hr_ovaries_January2015.fastq cleaned_ncRNA_matching_reads_TM3hs20hr_ovaries_January2015.txt).\n";
}
else
{
    $n=3;
    open(FASTQ,$ARGV[0]);
    while(<FASTQ>)
    {
	++$n;
	chomp;
	if ($n == 4)
	{
	    s/^@//;
	    s/ .*//;
	    $name=$_;
	    $n=0;
	}
        if ($n == 1)
        {
            $l1{$name}=$_;
        }
        if ($n == 3)
        {
            $l3{$name}=$_;
        }
    }
    close(FASTQ);

    open(LIST,$ARGV[1]);
    while (<LIST>)
    {
	chomp;
	print "\@$_\n$l1{$_}\n+\n$l3{$_}\n";
    }
    close(LIST);
}
