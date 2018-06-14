#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script arguments (developmental stage, and library number; e.g., ./Module_excludes_hairpin-matching_reads.pl embryon_36h 2).\n";
}
else
{
    open(BLACKSHEEP,"/mnt/data/home/herve.seitz/Amphioxus/miRNA_search/Hairpin-matching_reads_".$ARGV[0]."_".$ARGV[1].".sam");
    while(<BLACKSHEEP>)
    {
        if (/^\@PG\t/)
        {
            $read=1;
        }
        elsif ($read)
        {
            chomp;
            @array=split('\t',$_);
            if ($array[1]==0)
            {
                push(@suspect,$array[0]);
            }
        }
    }
    close(BLACKSHEEP);

    for $ori ('sense','antisense')
    {
        open(FA,$ori."_transcriptome-matching_reads_".$ARGV[0]."_".$ARGV[1].".fa");
        open(OUT,">$ori"."_transcriptome-matching_reads_not_matching_pre-miRNA_hairpins_".$ARGV[0]."_".$ARGV[1].".fa");
        while (<FA>)
        {
            chomp;
            if (/^>/)
            {
                s/^>//;
                $name=$_;
                if (grep {$_ eq $name} @suspect)
                {
                    $include=0;
                }
                else
                {
                    $include=1;
                }
            }
            else
            {
                if ($include)
                {
                    print OUT ">$name\n$_\n";
                }
            }
        }
    }
    close(OUT);
    close(FA);
}
