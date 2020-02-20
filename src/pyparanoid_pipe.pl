#!/bin/env perl

use strict;
use warnings;

my $threads = 78;
my $d = 'pan_flavo';

print "$d\n";
next if(-d "$d/gdb");
chdir $d;
mkdir "gdb" unless(-d 'gdb');
mkdir 'gdb/pep' unless(-d 'gdb/pep');
mkdir 'faa' unless(-d 'faa');
foreach my $fa (glob("*.fna")){
    my $pref = '';
    if($fa =~ m/(.+)\.f(n)?a(sta)?$/){
	$pref = $1;
	if($pref =~ m/\//){
	    my @parr = split(/\//, $pref);
	    $pref = $parr[(scalar(@parr)-1)];
	}
    }else{
	die "\"$fa\" does not have proper extension (.fasta, .fna, or .fa)\n";
    }
    unless(-e 'faa/'.$pref.'.faa'){
	system("prodigal -t $pref.tmp.ptrain -c -i $fa");
	system("prodigal -c -i $fa -a $pref.faa -d $pref.tmp.orf -t $pref.tmp.ptrain -f gff -o $pref.tmp.gff");
	system("mv *.faa faa");
	system("rm $pref.tmp*");
    }
}
open my $sfh, '>', 'strainlist.txt' or die $!;
foreach my $i (glob("faa/*.faa")){
    my @p = split(/\//, $i);
    my ($pref, $opref) = ($p[-1], $p[-1]);
    $opref =~ s/\.faa$/\.fna/;
    $pref =~ s/\.faa$//;
    my $o = 'gdb/pep/'.$pref.'.pep.fa';
    system("mv $i $o");
    system("cp ./$opref gdb/$pref.fasta");
    print $sfh "$pref\n";
}
close $sfh;
system("rm -r faa *.fna"); #OPTIONAL CLEANING STEP
system("BuildGroups.py --clean --verbose --use_MP --cpus $threads ./gdb strainlist.txt ./$d.pyp");
system("cp $d.pyp/strainlist.txt $d.pyp/prop_strainlist.txt");
system("cp $d.pyp/homolog.faa $d.pyp/prop_homolog.faa");
system("IdentifyOrthologs.py --use_MP --threshold 0.9 $d.pyp $d.ortho");

system("cp -r /home/mchevrette/jolab001/backup_workspace/20190730-bacillus_pyparanoid/pan_bacillus/mibig_gdb ./");
system("cp /home/mchevrette/jolab001/backup_workspace/20190730-bacillus_pyparanoid/pan_bacillus/mibig_strainlist.txt ./");
system("PropagateGroups.py mibig_gdb mibig_strainlist.txt ./$d.pyp");
chdir '..';
