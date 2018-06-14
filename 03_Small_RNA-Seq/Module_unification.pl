#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter name of input files (file containing blast output, then fata file that was used as an input for blast) (e.g., ./Module_unification.pl blast_output_unification.txt tmp_unified.fa).\n";
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
            s/ .*//;
            $name=$_;
        }   
        else
        {
            $seq{$name}=$seq{$name}.$_;
        }
    }
    close(FA);

    foreach $name (keys %seq)
    {
        $length{$name}=length($seq{$name});
    }

    open(BLAST,$ARGV[0]);
    while(<BLAST>)
    {
        chomp;
        @array=split(' ',$_);
        if (($array[0] ne $array[1]) && ($array[2]==100) && (($array[3]==length($seq{$array[0]})) || ($array[3]==length($seq{$array[1]}))))
        {
            $l0=length($seq{$array[0]});
            $l1=length($seq{$array[1]});
            if ($l0<=$l1)
            {
                $long=1;
                $short=0;
            }
            else
            {
                $long=0;
                $short=1;
            }
            $already=0;
            if (grep {$_ eq $array[$short]} @synonyms{$array[$long]})
            {
                $already=1;
            }
            if ($already==0)
            {
                push(@{$synonyms[$array[$long]]},$array[$short]);
                push(@{$synonyms[$array[$long]]},@{$synonyms[$array[$short]]});
                $merged{$array[$short]}=1;
            }
        }
    }
    close(BLAST);

    foreach $name (keys %seq)
    {
        if ($merged{$name} eq '')
        {
            print ">$name ; synonyms: @{$synonyms[$name]}\n$seq{$name}\n";    
        }
    }
} 
