#!/usr/bin/perl
if ($ARGV[0] eq '')
{
    print "Please enter input file name (e.g., ./Module_converts_fastq_to_fasta.pl 130404_SN132_A_L001_GZM-1_R1.fastq).\n";
}
else
{
    $file_name=$ARGV[0];
    $file_name=~s/\.fastq$//;
   
    open(OUT,">$file_name".".fa");
    open(IN,$ARGV[0]);
    while (<IN>)
    {
	++$n;
	if ($n == 1)
	{
	    s/^@/> /;
	    print OUT $_;
	}
	if ($n == 2)
	{
	    print OUT $_;
	}
	if ($n == 4)
	{
	    $n=0;
	}
    }
    close(IN);
    close(OUT);

    print "Done. Output file is named: $file_name.fa\n";
}
