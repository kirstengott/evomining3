#!/bin/env python

import sys, os, glob

## Usage: pyparanoid_pipeline.py <genome_directory>

## Set parameters
threads = 78
pyp_dir = 'pyp';
genome_dir = sys.argv[1]

## Set up pyparanoid input directories
if not os.path.exists(pyp_dir):
    os.mkdir(pyp_dir)
if not os.path.exists(pyp_dir+'/gdb'):
    os.mkdir(pyp_dir+'/gdb')
if not os.path.exists(pyp_dir+'/gdb/pep'):
    os.mkdir(pyp_dir+'/gdb/pep')

## Populate genomes list from the genome directory
genomes = []
for f in glob.glob(genome_dir+"/*.fasta"):
    genomes.append(f)
for f in glob.glob(genome_dir+"/*.fna"):
    genomes.append(f)
for f in glob.glob(genome_dir+"/*.fa"):
    genomes.append(f)
    
with open(pyp_dir+'/strainlist.txt', 'w') as sfh:
    for fna in genomes:    
        prefix = os.path.splitext(os.path.basename(fna))[0]
        ## Call genes with prodigal
        os.system(' '.join(['prodigal', '-t', prefix+'.tmp.ptrain', '-c', '-i', fna, '>', '/dev/null', '2>&1']))
        os.system(' '.join(['prodigal', '-c', '-i', fna, '-a', pyp_dir+'/gdb/pep/'+prefix+'.pep.fa', '-t', prefix+'.tmp.ptrain', '>', '/dev/null', '2>&1']))
        os.system(' '.join(['cp', fna, pyp_dir+'/gdb/'+prefix+'.fasta']))
        os.system(' '.join(['rm', prefix+'.tmp.ptrain']))
        sfh.write(prefix+"\n")

## Run pyparanoid
os.system(' '.join(['BuildGroups.py', '--clean', '--verbose', '--use_MP', '--cpus', str(threads), pyp_dir+'/gdb', pyp_dir+'/strainlist.txt', pyp_dir+'/pyparanoid']))
os.system(' '.join(['cp', pyp_dir+'/pyparanoid/strainlist.txt', pyp_dir+'/pyparanoid/prop_strainlist.txt']))
os.system(' '.join(['cp', pyp_dir+'/pyparanoid/homolog.faa', pyp_dir+'/pyparanoid/prop_homolog.faa']))
os.system(' '.join(['IdentifyOrthologs.py', '--use_MP', '--threshold', '0.9', pyp_dir+'/pyparanoid', pyp_dir+'/ortho']))

## Propogate to MIBiG
mibig_dir = os.path.realpath(__file__).replace('src/pyparanoid_pipeline.py', 'data/mibig_2.0')
mibig_list = os.path.realpath(__file__).replace('src/pyparanoid_pipeline.py', 'data/mibig_2.0_strainlist.txt')
os.system(' '.join(['PropagateGroups.py', mibig_dir, mibig_list, pyp_dir+'/pyparanoid']))
