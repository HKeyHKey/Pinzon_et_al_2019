#!/usr/bin/perl

$MIN_SIZE=18; # minimal read length
$MAX_SIZE=30; # maximal read length

if ($ARGV[0] eq '')
{
    print "Please enter script argument (e.g., ./Module_selects_from_blat_results.pl blat_output.psl).\n";
}
else
{
    open(IN,$ARGV[0]);
    while (<IN>)
    {
        chomp;
        @array=split('\t',$_);
        push(@{$info{$array[9]}},$array[19]);
        push(@{$block_nb{$array[9]}},$array[6]);
    }
    close(IN);

    foreach $transcript (keys %block_nb)
    {
        %seen=();
        @unique = do { %seen; grep { !$seen{$_}++ } @{$info{$transcript}} };
        $nb=@unique;
        if ($nb==1) #Selects transcripts whose exon/intron structure is the same among their genomic loci (if they map uniquely, no problem, but if they map at several places: the position of exon boundaries within the mRNA has to be the exact same)
        {
            if ($block_nb{$transcript}[0]>0) #Selects transcripts which have at least 1 intron
            {
                push(@selected_transcripts,$transcript);
                @{$block_starts{$transcript}}=split(',',$info{$transcript}[0]);
            }
        }
    }

    $nb_selected=@selected_transcripts;
    print "I will consider $nb_selected transcripts...\n";

    open(DETAILS,">Antisense_exon-exon_junction_read_details.dat");
    print DETAILS "Transcript Stage Library Read_ID Read_Sequence\n";

    foreach $stage ('embryon_8h','embryon_15h','embryon_36h','embryon_60h','male','femelle')
    {
        foreach $lib (1,2,3,4)
        {
            open(SAM,"Transcriptome_mapping_of_extragenomic_reads_".$stage."_".$lib.".sam");
            $read=0;
            while(<SAM>)
            {
                if (/^\@PG\t/)
                {
                    $read=1;
                }
                else
                {
                    if ($read)
                    {
                        chomp;
                        @array=split('\t',$_);
                        if (grep {$_ eq $array[2]} keys %block_starts)
                        {
                            $l=length($array[9]);
                            if (($l>=$MIN_SIZE) && ($l<=$MAX_SIZE))
                            {
                                $start=$array[3]-1; # Because SAM positions are 1-based and need to be converted to 0-based
                                $end=$start+$l-1;
                                $annot='exonic';
                                for $exon_start (@{$block_starts{$array[2]}})
                                {
                                    if (($start<$exon_start) && ($end>$exon_start))
                                    {
                                        $annot='junction';
### Below: if you find antisense exon-exon junction reads, report them in a separate file
                                        if ($array[1]==16)
                                        {
                                            $read_seq=$array[9];
                                            $read_seq=~tr/ACGT/TGCA/;
                                            $read_seq=reverse $read_seq;
                                            print DETAILS "$array[2] $stage $lib $array[0] $read_seq\n";
                                        }
### Above: if you find antisense exon-exon junction reads, report them in a separate file                                    
                                    }
                                }
                                ++$hit_count{$array[2]}{$annot}{$array[1]};
                            }
                        }
                    }
                }
            }
            close(SAM);
        }
    }
    close(DETAILS);
    
    $radical=$ARGV[0];
    $radical=~s/\.psl$//;
    open(OUT,">Junction_and_exonic_reads_$radical".".dat");
    print OUT "Transcript Exonic_reads_sense Junction_reads_sense Exonic_reads_antisense Junction_reads_antisense\n";
    for $transcript (keys %block_starts)
    {
        print OUT "$transcript";
        for $sense ('0','16')
        {
            for $annot ('exonic','junction')
            {
                if ($hit_count{$transcript}{$annot}{$sense})
                {
                    $display=$hit_count{$transcript}{$annot}{$sense};
                }
                else
                {
                    $display=0;
                }
                print OUT " $display";
            }
        }
        print OUT "\n";
    }
}
 
