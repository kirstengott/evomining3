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
focal_prefix = os.path.splitext(os.path.basename(focal_genome))[0]
focal_kofam = focal_prefix+'.pep.kofam_parsed.tsv'

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
ot = pd.read_csv('ortho_tall.tsv', delimiter="\t")

## Left join with in_groups and subset based on ANI criteria
ot = pd.merge(ot, i, 'left', on='strain_id')
ot = ot[ot.strain_id.isin(ani.strain_id)]

## Get functional data
func = pd.read_csv(focal_kofam, delimiter="\t")
func = func[['Gene', 'KO-D']]
gl = pd.read_csv('group_list.tsv', delimiter="\t")
gl = gl[['Gene', 'ortho_group']]
func = pd.merge(func, gl, 'left', on='Gene')
kegg_hier = os.path.realpath(__file__).replace('src/identify_expansions.py', 'data/kegg_hier.tsv')
kh = pd.read_csv(kegg_hier, delimiter="\t")
func = pd.merge(func, kh, 'left', on='KO-D')

## Check if ortho_groups are in metabolism and/or MIBiG
met_a = ['Metabolism'] ## Kegg level A matches for metabolism
vitaa_b = ['Amino acid metabolism', 'Metabolism of other amino acids', 'Metabolism of cofactors and vitamins'] ## Kegg level B matches for amino acid or vitamin metabolism
mb = pd.read_csv('mibig_groups.tsv', delimiter="\t")        
is_met = {}
is_vitaa = {}
is_mibig = {}
for index, row in func.iterrows():
    if row['DescA'] in met_a:
        is_met[row['ortho_group']] = 1
    if row['DescB'] in vitaa_b:
        is_vitaa[row['ortho_group']] = 1
    if row['ortho_group'] in mb.ortho_group.tolist():
        is_mibig[row['ortho_group']] = 1

## List of the group cutoffs
ingroups = list(i)[3:]

## Define the expansions
expansion = {}
hist = {}
## Loop through each cutoff group
for grp in ingroups:
    grpmed = {}
    ot_grp = ot.copy()
    ot_grp = ot_grp[[grp, 'ortho_group', 'n_ortho']]
    ot_grp[grp] = ot_grp[grp].astype(int)
    ## Summarize by ortho_group and cutoff group
    bygrp = ot_grp.groupby(['ortho_group', grp]).agg(['count', 'mean', 'std', 'median'])
    for index, row in bygrp.iterrows():
        ## Index is a tuple of (ortho_group, binary) where binary is 0 for not in cutoff group, 1 for in cutoff group
        if index[0] not in grpmed:
            grpmed[index[0]] = {}
        ## In cutoff group
        if index[1]==1:
            grpmed[index[0]]['Yes'] = row[('n_ortho', 'median')]
        ## Not in cutoff group
        else:
            grpmed[index[0]]['No'] = row[('n_ortho', 'median')]
    ## Save to dictionaries
    for o in grpmed:
        if not o in expansion:
            expansion[o] = {}
        expansion[o][grp] = grpmed[o]['Yes'] - grpmed[o]['No']
        if not grp in hist:
            hist[grp] = []
        hist[grp].append(grpmed[o]['Yes'] - grpmed[o]['No'])

## Save data externally
with open('expansions_list.tsv', 'w') as efh:
    efh.write("\t".join(['ortho_group']+ingroups+['is_metabolism', 'is_vit_or_aa_metabolism', 'is_mibig'])+"\n")
    for ortho_group in sorted(expansion):
        efh.write(ortho_group)
        for grp in ingroups:
            efh.write("\t"+str(expansion[ortho_group][grp]))
        efh.write("\t")
        if ortho_group in is_met:
            efh.write('Yes')
        else:
            efh.write('No')
        efh.write("\t")
        if ortho_group in is_vitaa:
            efh.write('Yes')
        else:
            efh.write('No')
        efh.write("\t")
        if ortho_group in is_mibig:
            efh.write('Yes')
        else:
            efh.write('No')
        efh.write("\n")

## Generate histograms for each cutoff group
for grp in ingroups:
    hist_plot = sb.distplot(hist[grp], bins=np.arange(-20,20.5,.5), kde=False)
    png = hist_plot.get_figure()
    png.axes[0].set_yscale('log')
    png.savefig(grp+'.hist.png')
    plt.close()

