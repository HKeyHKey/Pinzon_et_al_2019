#!/usr/bin/perl

$MIN_LENGTH=100; #minimal length (in codons) for an ORF to be considered
$MAX_RANK=3; #maximal rank of the AUG to be considered (here: don't scan past the 3rd AUG in the transcript)
$MIN_POSITION=1; #minimal coordinate of the AUG in the transcript (0-based; here: exclude transcripts with 0 nt long 5´ UTRs)

if ($ARGV[0] eq '')
{
    print "Please give script argument (e.g., ./Module_ORF_detector.pl Bla_annot_final_refTranscripts.fa).\n";
}
else
{
    open(DATA,$ARGV[0]);
    while(<DATA>)
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
    close(DATA);

    open(OUT1,">ORF_length_data.dat");
    print OUT1 "Transcript_name Start_codon_position Start_codon_rank ORF_length\n";
    $file_name=$ARGV[0];
    $file_name=~s/\.fa$//;
    open(SELECT,">Transcripts_with_convicing_ORFs_in_$file_name".".txt");
    open(UNSELECT,">Transcripts_without_convicing_ORFs_in_$file_name".".txt");
    foreach $name (keys %seq)
    {
#	print "Now analyzing $name...\n";
	$select=0;
	$offset=0;
	$start_index=0;
        $start=index($seq{$name},'ATG',$offset);
#	print "   start=$start\n";
        while ($start != -1)
        {
	    ++$start_index;
	    $nt=$start;
	    $stop='';
	    while (($nt+2<length($seq{$name})) && ($stop eq ''))
	    {
		$extract=substr $seq{$name},$nt,3;
		if (($extract eq 'TAA') || ($extract eq 'TAG') || ($extract eq 'TGA'))
		{
		    $stop=$nt;
		}
		$nt=$nt+3;
	    }
	    if ($stop eq '')
	    {
		$stop=$nt;#if I could not find a stop codon, I consider that it is immediately after the end of the sequence
	    }
	    $ORF_length=($stop-$start)/3;
	    if (($ORF_length>=$MIN_LENGTH) && ($start_index<=$MAX_RANK) && ($start>=$MIN_POSITION))
	    {
		$select=1;
	    }
	    print OUT1 "$name $start $start_index $ORF_length\n";
            $offset=$start+1;
	    $start=index($seq{$name},'ATG',$offset);
        }
	if ($select)
	{
	    print SELECT "$name\n";
	}
	else
	{
	    print UNSELECT "$name\n";
	}
    }
    close(SELECT);
    close(UNSELECT);
}
