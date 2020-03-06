#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

## USAGE STATEMENT
if (!defined($ARGV[0]) or $ARGV[0] eq '-help') {
die ("
Program: as5 helper

Usage:    ./as500_pipe.pl [options] fasta_in_files

Options:
     -as5_exec            path to antismash5 docker executable
     -taxon               {bacteria, fungi}
     -genefindingtool     {glimmerhmm,prodigal,prodigal-m,none,error}
     -threads             number of threads
     -help                Print this usage statement

");
}

my $taxon = 'bacteria';
my $genefindingtool = 'prodigal';
my $threads = 80;
my $cmd = 'run_antismash500';

GetOptions(
    'as5_exec' => \$cmd,
    'taxon' => \$taxon,
    'genefindingtool' => \$genefindingtool,
    'threads' => \$threads
    );


#my $args = '--genefinding-tool prodigal --cpus '.$threads.' --minimal'; ## Minimal
my $args = '--taxon '.$taxon.' --genefinding-tool '.$genefindingtool.' --cpus '.$threads.' --clusterhmmer --smcog-trees --cb-general --cb-subclusters --cb-knownclusters --asf --pfam2go --cb-nclusters 50'; ## Do it all
my $d = 'as500';
mkdir $d unless(-d $d);
my $tot = scalar(@ARGV);
my $n = 1;
foreach my $fa (@ARGV){
    my $p = '';
    if($fa =~ m/(.+)\.f(n)?a(sta)?$/){
	$p = $1;
	if($p =~ m/\//){
	    my @parr = split(/\//, $p);
	    $p = $parr[-1];
	}
    }else{
	die "\"$fa\" does not have proper extension (.fasta, .fna, or .fa)\n";
    }
    ## Run antismash
    my $asd = './'.$d.'/' . $p . '_as500';
    #system("rm -r $asd") if(-d $asd);
    unless(-d $asd){
	print "Running $p antiSMASH 5.0.0 ($n of $tot)...";
	system("$cmd $fa $asd $args");
    }
    print "DONE!\n";
    $n += 1;
}
