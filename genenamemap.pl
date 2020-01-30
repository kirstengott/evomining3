#!/bin/env perl

use strict;
use warnings;


open my $gfh, '<', 'pan_flavo/pan_flavo.pyp/homolog.faa' or die $!;
print join("\t", 'genome', 'gene_name', 'contig', 'contig_gene_num', 'ortho_group')."\n";
while(<$gfh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	my($g, $gn, $grp) = split(/\|/, $1);
	my $contig = $gn;
	$contig =~ s/_\d+$//;
	my $cgnum = $gn;
	$cgnum =~ s/.+_(\d+)$/$1/;
	print join("\t", $g, $gn, $contig, $cgnum, $grp)."\n";
    }
}
close $gfh;
