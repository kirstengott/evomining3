#!/bin/env perl

use strict;
use warnings;

## Usage: perl calculate_ani.pl focalgenome.fna genome_folder

## Set internal parameters
my $threads = 60;

## Parse command line paramaters
my ($focal_genome, $genome_directory) = (shift, shift);

## Check focal genome is correct file extension
die "\"$focal_genome\" does not have proper extension (.fasta, .fna, or .fa)\n" unless($focal_genome =~ m/(.+)\.f(n)?a(sta)?$/);
my @focal_path = split(/\//, $focal_genome);
my $focal_prefix = $focal_path[-1];
$focal_prefix =~ s/\.f(n)?a(sta)?$//;

## Open output file handle and write header
open my $ofh, '>', 'ANI.tsv' or die $!;
print $ofh join("\t", 'Query', 'Reference', 'ANI', 'Mappings', 'TotalFrags', 'QueryAligned')."\n";

## Loop through each query genome in the genome directory
foreach my $query (glob("$genome_directory/*.fna"), glob("$genome_directory/*.fasta"), glob("$genome_directory/*.fa")){
    next if($query eq $focal_genome);
    my @query_path = split(/\//, $query);
    my $query_prefix = $query_path[-1];
    $query_prefix =~ s/\.f(n)?a(sta)?$//;
    ## Perform ANI calculation
    system("fastANI -q $query -r $focal_genome --visualize -t $threads -o ani_tmp.txt > /dev/null 2>&1");
    ## Parse results to get the length of ANI alignment
    open my $vfh, '<', 'ani_tmp.txt.visual' or die $!;
    my $len = 0;
    while(<$vfh>){
	chomp;
	my @ln = split(/\t/, $_);
	$len += abs($ln[6] - $ln[7])+1;
    }
    close $vfh;
    ## Parse results to get the ANI summary information
    my ($ani, $map, $tf) = (undef, undef, undef);
    open my $tfh, '<', 'ani_tmp.txt' or die $!;
    while(<$tfh>){
	chomp;
	my @ln = split(/\t/, $_);
	($ani, $map, $tf) = ($ln[2], $ln[3], $ln[4]);
    }
    close $tfh;
    ## If nothing, then set to zero
    $ani = 0 unless(defined($ani) && $ani =~ m/\d+/);
    $map = 0 unless(defined($map) && $map =~ m/\d+/);
    $tf = 0 unless(defined($tf) && $tf =~ m/\d+/);
    print $ofh join("\t", $query_prefix, $focal_prefix, $ani, $map, $tf, $len)."\n";
    ## Clean up
    system("rm ani_tmp.txt*");
}
close $ofh;							 
