#!/bin/env python

import glob, os, re
from Bio import SeqIO

as5 = {}
with open('as500.tsv', 'r') as afh:
    for ln in afh.read().splitlines():
        genome, contig, start, end = ln.split("\t")
        if genome not in as5:
            as5[genome] = {}
        if contig not in as5[genome]:
            as5[genome][contig] = {}
        as5[genome][contig][int(start)] = int(end)

print("\t".join(['genome', 'contig', 'gene_id', 'start', 'end', 'inclust']))
for faa in glob.glob("pan_flavo/gdb/pep/*.pep.fa"):
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
        print("\t".join([genome, contig, seq.id, header[1], header[2], inclust]))
        i += 1
