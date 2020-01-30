#!/bin/env perl

use strict;
use warnings;

print join("\t", 'strain_id', 'Gene', 'ortho_group')."\n";
open my $hfh, '<', 'pan_flavo/pan_flavo.pyp/homolog.faa' or die $!;
while(<$hfh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	my @l = split(/\|/, $1);
	print join("\t", @l)."\n" if($l[0] == '89329');
    }
}
close $hfh;
