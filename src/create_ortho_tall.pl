#!/bin/env perl

use strict;
use warnings;

my $hmat = 'pan_flavo/pan_flavo.pyp/homolog_matrix.txt';
my @header = ();
my %mibig = ();
open my $hfh, '<', $hmat or die $!;
open my $ofh, '>', 'ortho_tall.tsv' or die $!;
print $ofh join("\t", 'ortho_group', 'strain_id', 'n_ortho')."\n";
while(<$hfh>){
    chomp;
    if(scalar(@header)==0){
	$_ =~ s/^\s+//;
	@header = split(/\t/, $_);
	next;
    }else{
	my ($group, @ln) = split(/\t/, $_);
	my %seen = ();
	for(my $i=0;$i<scalar(@header);$i+=1){
	    if($header[$i] =~ m/^BGC/ && $ln[$i]>0){
		$mibig{$group} += 1; 
	    }else{
		print $ofh join("\t", $group, $header[$i], $ln[$i]/2)."\n" unless(exists $seen{$header[$i]});
		$seen{$header[$i]} += 1;
	    }
	}
    }
}
close $hfh;
close $ofh;
open my $mfh, '>', 'mibig_groups.tsv' or die $!;
print $mfh join("\t", 'ortho_group', 'count')."\n";
foreach my $group (keys %mibig){
    print $mfh join("\t", $group, $mibig{$group})."\n"
}
close $mfh;
