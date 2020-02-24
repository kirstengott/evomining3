#!/bin/env python

import os, sys, re

## Usage: python kofam_annotation <ko_list> <prokaryote.hal> <focal_genome_pep.faa>

## Note: this is part of kofamscan and a *hal database must be downloaded from
## ftp://ftp.genome.jp/pub/tools/kofamscan/
## with the corresponding *hmms
ko_list, hal_db = sys.argv[1], sys.argv[2]
focal_pep = sys.argv[3]

## Functional annotation with kofamscan
threads = 12
focal_prefix = os.path.splitext(os.path.basename(focal_pep))[0]
kofam_out = focal_prefix+'.kofam.txt'
if not os.path.isfile(kofam_out):
    os.system(' '.join(['exec_annotation', '-k', ko_list, '-p', hal_db, '--cpu='+str(threads), '-o', kofam_out, focal_pep]))

## Parse the results
kofam_parsed = focal_prefix+'.kofam_parsed.tsv'
with open(kofam_parsed, 'w') as kfh:
    kfh.write("\t".join(['Genome', 'Gene', 'KO-D'])+"\n")
    with open(kofam_out, 'r') as rfh:
            for ln in rfh.read().splitlines():
                if ln.startswith('#'):
                    continue
                if ln.startswith('*'):
                    ln = ln[1:]
                ln = ln.strip()
                ln = re.sub('\s+', ' ', ln)
                gene, ko, thresh, score, evalue, desc = ln.split(' ', 5)
                if thresh == '-':
                    thresh = str(0)
                if float(score) >= float(thresh):
                    kfh.write("\t".join([focal_prefix, gene, ko])+"\n")
