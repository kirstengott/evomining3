#!/bin/env python

import os, sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sb
import numpy as np

## Usage: python identify_ani_groups.py <focal_genome_fasta>

## Define alignment fraction:
## (minimum number of bp in alignment) / (bp in focal genome)
aln_frac = 0.25

## Parse command line paramaters
focal_genome = sys.argv[1]

## Get focal genome length and calculate alignment cutoff
focal_len = 0
with open(focal_genome, 'r') as ffh:
    for ln in ffh.read().splitlines():
        if not ln.startswith('>'):
            focal_len += len(ln)
aln_cut = aln_frac * focal_len

## Parse ANI file
ani = pd.read_csv('ANI.tsv', delimiter="\t")
ani = ani[ani['ANI'] > 0]
ani = ani[ani['QueryAligned'] > aln_cut]
ani = ani[['Query', 'ANI']]
ani.columns = ['strain_id', 'ANI']

## Parse group and ortho files
i = pd.read_csv('in_groups.tsv', delimiter="\t")
mb = pd.read_csv('mibig_groups.tsv', delimiter="\t")
ot = pd.read_csv('ortho_tall.tsv', delimiter="\t")

ot = pd.merge(ot, i, 'left', on='strain_id')
ot = ot[ot.strain_id.isin(ani.strain_id)]

ingroups = list(i)[3:]

expansion = {}
hist = {}
for grp in ingroups:
    grpmed = {}
    ot_grp = ot.copy()
    ot_grp = ot_grp[[grp, 'ortho_group', 'n_ortho']]
    ot_grp[grp] = ot_grp[grp].astype(int)
    bygrp = ot_grp.groupby(['ortho_group', grp]).agg(['count', 'mean', 'std', 'median'])
    for index, row in bygrp.iterrows():
        if index[0] not in grpmed:
            grpmed[index[0]] = {}
        if index[1]==1:
            grpmed[index[0]]['Yes'] = row[('n_ortho', 'median')]
        else:
            grpmed[index[0]]['No'] = row[('n_ortho', 'median')]
    for o in grpmed:
        if not o in expansion:
            expansion[o] = {}
        expansion[o][grp] = grpmed[o]['Yes'] - grpmed[o]['No']
        if not grp in hist:
            hist[grp] = []
        hist[grp].append(grpmed[o]['Yes'] - grpmed[o]['No'])

with open('expansions_list.tsv', 'w') as efh:
    efh.write("\t".join(['ortho_group']+ingroups)+"\n")
    for ortho_group in sorted(expansion):
        efh.write(ortho_group)
        for grp in ingroups:
            efh.write("\t"+str(expansion[ortho_group][grp]))
        efh.write("\n")

for grp in ingroups:
    hist_plot = sb.distplot(hist[grp], bins=np.arange(-20,20.5,.5), kde=False)
    png = hist_plot.get_figure()
    png.savefig(grp+'.hist.png')
    plt.close()
