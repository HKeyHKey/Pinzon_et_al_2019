#!/usr/bin/perl

$GENOME='/mnt/data/home/herve.seitz/Genomes/B_lanceolatum/assembly/Bl71nemr.fa';
$SIGNIF_CUTOFF=0.000001;

if ($ARGV[0] eq '')
{
    print "Please enter name of file containing blast output (e.g., ./Module_extracts_miRNA_hairpins.pl blast_output_hairpins.txt).\n";
}
else
{
    open(GENOME,$GENOME);
    while(<GENOME>)
    {
        chomp;
        if (/^>/)
        {
            s/^> *//;
            s/ .*//;
            $name=$_;
        }   
        else
        {
            $seq{$name}=$seq{$name}.$_;
        }
    }
    close(GENOME);

    open(OUT,">Orthologs_to_known_hairpins.fa");
    open(IN,$ARGV[0]);
    while(<IN>)
    {
        chomp;
        @array=split('\t',$_);
        if ($array[10]<=$SIGNIF_CUTOFF)
        {
            ($ortholog,$scaffold,$start,$end)=($array[0],$array[1],$array[8],$array[9]);
            if ($start<$end)
            {
                $extract=substr $seq{$scaffold},$start,$end-$start+1;
                $strand='+';
            }
            else
            {
                $extract=substr $seq{$scaffold},$end,$start-$end+1;
                $extract=reverse $extract;
                $extract=~tr/ACGTacgt/TGCAtgca/;
                $strand='-';
            }
            print OUT ">$ortholog ortholog at $scaffold bp $start-$end ($strand strand)\n$extract\n";
        }
    }
    close(IN);
    close(OUT);
}   
