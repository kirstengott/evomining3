#!/bin/env python

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sb

## Usage: python identify_ani_groups.py <ANIfile>

ani_file = 'ANI.tsv'

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
ani = pd.read_csv(ani_file, delimiter="\t")
ani = ani[ani['ANI'] > 0]
ani = ani[ani['QueryAligned'] > aln_cut]

## Plot the histogram
hist_plot = sb.distplot(ani['ANI'], bins=40)
png = hist_plot.get_figure()
png.savefig("ani_histogram.png")
