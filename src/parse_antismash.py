#!/bin/env python

import os, re, sys, glob
from Bio import SeqIO

## Usage: python parse_antismash.py <pyparanoid_pep_folder> <antismash_final_gbk_1> <antismash_final_gbk_2>....<antismash_final_gbk_n>

## Read in antiSMASH BGC genome coordinates
as5 = {}
with open('antismash.tsv', 'w') as afh:
    for gbk in sys.argv[2:]:
        if not re.search( r"region\d+", gbk ): ## We only want the genome-wide gbk, so ignore orthers
            genome = os.path.split(gbk)[0].split('/')[-1]
            for seq in SeqIO.parse(gbk, "genbank"):
                for feat in seq.features:
                    if feat.type == 'region':
                        afh.write("\t".join([genome, seq.id, str(feat.location.start), str(feat.location.end)])+"\n")
                        if genome not in as5:
                            as5[genome] = {}
                        if seq.id not in as5[genome]:
                            as5[genome][seq.id] = {}
                        as5[genome][seq.id][feat.location.start] = feat.location.end

## Determine which genes are within antiSMASH-defined BGCs
with open('gene_in_bgc.tsv', 'w') as gfh:
    gfh.write("\t".join(['genome', 'contig', 'gene_id', 'start', 'end', 'inclust'])+"\n")
    for faa in glob.glob(sys.argv[1]+"/*.pep.fa"):
        genome = os.path.split(faa)[-1].split('.')[0]
        i = 0
        for seq in SeqIO.parse(faa, "fasta"):
            header = re.split('\s#\s', seq.description)
            contig = '_'.join(seq.id.split('_')[0:-1])
            inclust = 'False'
            if genome in as5:
                if contig in as5[genome]:
                    for as_start in as5[genome][contig]:
                        if int(header[1]) >= as_start and int(header[1]) <= as5[genome][contig][as_start]:
                            inclust = 'True'
                        elif int(header[2]) >= as_start and int(header[2]) <= as5[genome][contig][as_start]:
                            inclust = 'True'
            gfh.write("\t".join([genome, contig, seq.id, header[1], header[2], inclust])+"\n")
            i += 1
