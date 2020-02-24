#!/bin/env python

import sys, os

## Usage: python create_ortho_tall.py <focal_genome.fasta>

## Define inputs
pyp_dir = 'pyp/pyparanoid'
homolog_faa = pyp_dir+'/homolog.faa'
homolog_mat = pyp_dir+'/homolog_matrix.txt'
focal_prefix = os.path.splitext(os.path.basename(sys.argv[1]))[0]

## Generate the ortholog group list for the focal genome
with open('group_list.tsv', 'w') as gfh:
    gfh.write("\t".join(['strain_id', 'Gene', 'ortho_group'])+"\n")
    with open(homolog_faa, 'r') as ffh:
        for ln in ffh.read().splitlines():
            if ln.startswith('>'):
                sid, gene, og = ln[1:].split('|')
                if sid == focal_prefix:
                    gfh.write("\t".join([sid, gene, og])+"\n")

## Read in the homolog matrix and generate the tall ortholog table
header = []
mibig = {}
with open('ortho_tall.tsv', 'w') as ofh:
    ofh.write("\t".join(['ortho_group', 'strain_id', 'n_ortho'])+"\n")
    with open(homolog_mat, 'r') as mfh:
        for ln in mfh.read().splitlines():
            if len(header) == 0:
                header = ln[1:].split("\t")
            else:
                seen = {}
                group, counts = ln.split("\t", 1)
                for i, count in enumerate(counts.split("\t")):
                    if header[i].startswith('BGC'):
                        if float(count) > 0:
                            if group in mibig:
                                mibig[group] += 1
                            else:
                                mibig[group] = 1
                    else:
                        if header[i] not in seen:
                            ofh.write("\t".join([group, header[i], str(float(count)/2)])+"\n")
                            seen[header[i]] = 1

with open('mibig_groups.tsv', 'w') as mgfh:
    mgfh.write("\t".join(['ortho_group', 'n_ortho'])+"\n")
    for group in mibig:
        mgfh.write("\t".join([group, str(mibig[group])])+"\n")
