#!/usr/bin/perl

if ($ARGV[0] eq '')
{
    print "Please enter input file name (e.g., ./Module_parse_XML.pl XK9RHVPE01R-Alignment_Single_file.xml).\n";
}
else
{    
    open(IN,$ARGV[0]);
    while(<IN>)
    {
	chomp;
	if (/^ *<query-title>/)
	{
	    $query=$_;
	    $query=~s/^ *<query-title>//;
	    $query=~s/<\/query-title>$//;
	}
	if (/^ *<title>/)
	{
	    $subject=$_;
	    $subject=~s/^ *<title>//;
	    $subject=~s/<\/title>$//;
	}
	if (/^ *<evalue>/)
	{
	    $evalue=$_;
	    $evalue=~s/^ *<evalue>//;
	    $evalue=~s/<\/evalue>$//;
	    if ($e{$query}{$subject} eq '')
	    {
		$e{$query}{$subject}=$evalue; #Fill the hash with the first recorded E-value for that query/subject pair (that's the lowest E-value for that pair)
	    }
	}
    }
    close(IN);

    $radical=$ARGV[0];
    $radical=~s/\..*//;
    open(OUT,">Annotated_best_hits_from_$radical".".txt");
    print OUT "Contig\tBest hit\tE-value for its best HSP\n";
    open(DETAILED,">Detailed_annotated_best_hits_from_$radical".".txt");
    print DETAILED "Contig\tBest hit\tE-value for its best HSP\n";
    for $query (keys %e)
    {
	$min='';
	$best='';
	for $subject (keys %{$e{$query}})
	{
	    if (($min eq '') || ($e{$query}{$subject}<=$min))
	    {
		$min=$e{$query}{$subject};
		$best=$best.'/'.$subject;
	    }
	}
	$best=~s/^\///;
	print DETAILED "$query\t$best\t$min\n";

	@array=split('/',$best);
	@simplified=();
	foreach $element (@array)
	{
	    $element=~s/ /_/;
	    $element=~s/ .*//;
	    $element=~s/_/ /;
	    push(@simplified,$element);
	}
	%seen=();
	@unique = do {%seen; grep { !$seen{$_}++ } @simplified };
	$display=join('/',@unique);
	print OUT "$query\t$display\t$min\n";
    }
    close(OUT);
    close(DETAILED);
}
