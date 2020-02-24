#!/bin/env python

import sys, os, glob

## Usage: python calculate_ani.py focalgenome.fna genome_folder

## Set internal parameters
threads = 60;

## Parse command line paramaters
focal_genome, genome_dir = sys.argv[1], sys.argv[2]

## Check focal genome is correct file extension
focal_ext = os.path.splitext(focal_genome)[-1]
if focal_ext != '.fasta' and focal_ext != '.fna' and focal_ext != '.fa':
    raise Exception("focal_genome does not have proper extension (.fasta, .fna, or .fa)\n")
focal_prefix = os.path.splitext(os.path.basename(focal_genome))[0]

## Populate queries list from the genome directory
queries = []
for f in glob.glob(genome_dir+"/*.fasta"):
    queries.append(f)
for f in glob.glob(genome_dir+"/*.fna"):
    queries.append(f)
for f in glob.glob(genome_dir+"/*.fa"):
    queries.append(f)

## Open output file handle and write header
with open('ANI.tsv', 'w') as afh:
    afh.write("\t".join(['Query', 'Reference', 'ANI', 'Mappings', 'TotalFrags', 'QueryAligned'])+"\n")
    ## Loop through each query genome
    for query in queries:
        query_prefix = os.path.splitext(os.path.basename(query))[0]
        ## Perform ANI calculation
        os.system(' '.join(['fastANI', '-q', query, '-r', focal_genome, '--visualize', '-t', str(threads), '-o', 'ani_tmp.txt', '>', '/dev/null', '2>&1']))
        ## Parse results to get the length of ANI alignment
        alilen = 0
        with open('ani_tmp.txt.visual', 'r') as vfh:
            for ln in vfh.read().splitlines():
                col = ln.split("\t")
                alilen += abs(float(col[6]) - float(col[7]))+1;
        ## Parse results to get the ANI summary information
        ani, maps, tf = 0, 0, 0
        with open('ani_tmp.txt', 'r') as tfh:
            for ln in tfh.read().splitlines():
                col = ln.split("\t")
                ani, maps, tf = col[2], col[3], col[4]
        afh.write("\t".join([query_prefix, focal_prefix, str(ani), str(maps), str(tf), str(alilen)])+"\n")
        for tmp in ['ani_tmp.txt', 'ani_tmp.txt.visual']:
            os.system(' '.join(['rm', tmp]))
