#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please give script arguments (e.g., ./Module_TXT_to_FA.pl Transcripts_with_convicing_ORFs_in_Bla_annot_final_refTranscripts.txt Bla_annot_final_refTranscripts.fa).\n";
}
else
{
    open(FA,$ARGV[1]);
    while(<FA>)
    {
	chomp;
	if (/^>/)
	{
	    s/^> *//;
	    $name=$_;
	}
	else
	{
	    $seq{$name}=$seq{$name}.$_;
	}
    }
    close(FA);

    open(TXT,$ARGV[0]);
    while(<TXT>)
    {
        chomp;
        print ">$_\n$seq{$_}\n";
    }
    close(TXT);
}
