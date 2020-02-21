#!/usr/bin/perl

use strict;
use warnings;

#my $in_dir="../";
my $in_file=$ARGV[0];

#--------------------------------------------------
open my $sfh, '<', $in_file or die $!;
my @genomes = ();
while(<$sfh>){
    chomp;
    my @ln = split(/\t/, $_);
	print "$ln[0]\n";
    push @genomes, 'fasta/'.$ln[0].'.fasta';
}
close $sfh;
#---------------------------------------------------------------------------------
open my $ofh, '>', 'ANI_fj.tsv' or die $!;
print $ofh join("\t", 'Query', 'Reference', 'Genome1', 'Genome2', 'ANI', 'Mappings', 'TotalFrags', 'QueryAligned')."\n";
foreach my $q (@genomes){
    my $qpref = $q;
    $qpref =~ s/.+\/(.+)\.fasta$/$1/;
    next unless($qpref eq '133318');
    foreach my $ref (@genomes){
        next if($q eq $ref);
	my $rpref = $ref;
	$rpref =~ s/.+\/(.+)\.fasta$/$1/;
        system("fastANI -q $q -r $ref --visualize -t 60 -o fjtmp.txt");
	open my $vfh, '<', 'fjtmp.txt.visual' or die $!;
	my $len = 0;
        while(<$vfh>){
            chomp;
            my @ln = split(/\t/, $_);
            $len += abs($ln[6] - $ln[7])+1;
        }
        close $vfh;
	my ($ani, $map, $tf) = (undef, undef, undef);
	open my $tfh, '<', 'fjtmp.txt' or die $!;
	while(<$tfh>){
	    chomp;
	    my @ln = split(/\t/, $_);
	    ($ani, $map, $tf) = ($ln[2], $ln[3], $ln[4]);
	}
	close $tfh;
	$ani = 0 unless($ani =~ m/\d+/);
	$map = 0 unless($map =~ m/\d+/);
	$tf = 0 unless($tf =~ m/\d+/);
	my @alph = sort ($qpref, $rpref);
	print $ofh join("\t", $qpref, $rpref, @alph, $ani, $map, $tf, $len)."\n";
        system("rm fjtmp.txt*");
    }
}
close $ofh;
