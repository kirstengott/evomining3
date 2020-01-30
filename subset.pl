#!/bin/env perl

use strict;
use warnings;

my %keep = ();
open my $ifh, '<', 'in_groups.tsv' or die $!;
while(<$ifh>){
    chomp;
    my ($r, @rest) = split(/\t/, $_);
    next if($r eq 'Reference');
    $keep{$r} = 1;
}
close $ifh;
foreach my $fna (glob("pan_flavo/*.fna")){
    my $pref = $fna;
    $pref =~ s/^pan_flavo\/(.+)\.fna$/$1/;
    system("rm $fna") unless(exists $keep{$pref});
}
