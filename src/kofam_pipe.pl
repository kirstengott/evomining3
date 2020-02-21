#!/bin/env perl

use strict;
use warnings;

my $cmd = '/home/mchevrette/builds/kofamscan/bin/exec_annotation';
foreach my $faa (@ARGV){
    my @p = split(/\//, $faa);
    my $pref = $p[-1];
    $pref =~ s/\.faa$//;
    my $txt = $pref.'.kofam.txt';
    system("$cmd -o $txt $faa");
    open my $ofh, '>', $pref.'.kofam_parsed.tsv' or die $!;
    print $ofh join("\t", 'Genome', 'Gene', 'KO-D')."\n";
    open my $tfh, '<', $txt or die $!;
    while(<$tfh>){
        chomp;
        next if($_ =~ m/^#/);
        $_ =~ s/^\**\s+//;
        my ($gene, $ko, $thresh, $score, $evalue, $desc) = split(/\s+/, $_, 6);
        $thresh = 0 if($thresh eq '-');
        print $ofh join("\t", $pref, $gene, $ko)."\n" if($score >= $thresh);
    }
    close $tfh;
    close $ofh;
}
