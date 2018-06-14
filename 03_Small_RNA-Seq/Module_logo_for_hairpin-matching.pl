#!/usr/bin/perl

if ($ARGV[0] eq '')
{
    print "Please enter name of input SAM file (e.g., ./Module_size_dist_for_hairpin-matching.pl Hairpin-matching_reads_embryon_36h_4.sam).\n";
}
else
{
    open(IN,$ARGV[0]);
    while(<IN>)
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
                $l=length $array[9];
                $seq=$array[9];
                $seq=~tr/T/U/;
                $sense_seq{$l}{$array[0]}=$seq;
            }
            if ($array[1]==16)
            {
                $l=length $array[9];
                $seq=$array[9];
                $seq=~tr/T/U/;
                $antisense_seq{$l}{$array[0]}=$seq;
            }
        }
    }
    close(IN);

    ($sample,$lib)=($ARGV[0],$ARGV[0]);
    $sample=~s/Hairpin-matching_reads_//;
    $sample=~s/_\d\.sam$//;
    $lib=~s/.*_//;
    $lib=~s/\.sam$//;
    for ($size=18;$size<=30;++$size)
    {
        $n{'sense'}=0;
        open(OUT,">$size"."-mers_sense_hairpin-matching_reads_$sample"."_$lib".".fa");
        for $id (keys %{$sense_seq{$size}})
        {
            print OUT ">$id\n$sense_seq{$size}{$id}\n";
            ++$n{'sense'};
        }
        close(OUT);
        $n{'antisense'}=0;
        open(OUT,">$size"."-mers_antisense_hairpin-matching_reads_$sample"."_$lib".".fa");
        for $id (keys %{$antisense_seq{$size}})
        {
            print OUT ">$id\n$antisense_seq{$size}{$id}\n";
            ++$n{'antisense'};
        }
        for $ori ('sense','antisense')
        {
            if ($n{$ori})
            {
                `/mnt/data/home/herve.seitz/Annotation_pipeline/weblogo/seqlogo -f $size"-mers_"$ori"_hairpin-matching_reads_"$sample"_"$lib".fa" -F EPS -o Logo_$size"-mers_"$ori"_hairpin-matching_reads_"$sample"_"$lib -abcMnY`;
            }
        }
    }
}   
