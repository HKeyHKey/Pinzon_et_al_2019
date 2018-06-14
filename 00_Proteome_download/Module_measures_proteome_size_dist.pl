#!/usr/bin/perl
if ($ARGV[0] eq '')
{
    print "Please give script arguments (e.g., ./Module_measures_proteome_size_dist.pl Octodon_degus_proteome.fa)."."\n";
}
else
{
    $file=$ARGV[0];
    open(HITS,$file);
    while (<HITS>)
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
    close(HITS);
    foreach $name (keys %seq)
    {
        $l=length($seq{$name});
        print "$name $l\n";
    }
}

