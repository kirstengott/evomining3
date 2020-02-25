#!/bin/env python

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sb
import numpy as np

## Usage: python identify_ani_groups.py

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

## Identify adjacent zero values in the histogram
## The average of these valleys is used to estimate group breakpoints
## Bin width is set to 0.1 and all ANI are rounded to one decimal
ani_hist = {}
for a in ani.ANI.tolist():
    onedec = round(a, 1)
    if onedec in ani_hist:
        ani_hist[onedec] += 1
    else:
        ani_hist[onedec] = 1
breaks = []
i = min(ani_hist.keys()) - 0.1
zero_flag = False ## Are we in a sequence of zeroes?
first_zero = i ## Where did we start? Initialized on i
len_zero = 0 ## How many zeroes?
min_zero = 3 ## Minimum allowed length of a sequence of zeroes to call a break
while i <= 100:
    i = round(i,1)
    if i in ani_hist:
        if zero_flag:
            if len_zero >= min_zero:
                last_zero = i - 0.1
                breaks.append(round( (first_zero+last_zero) / 2, 2))
            zero_flag = False
            len_zero = 0
    else:
        len_zero += 1
        if not zero_flag:
            first_zero = i
            zero_flag = True
    i += 0.1
print("Suggested break points to define ANI groups are:\t"+str(breaks))

## Plot the histogram
hist_plot = sb.distplot(ani['ANI'], bins=40)
pdf = hist_plot.get_figure()
for b in breaks:
    plt.axvline(x=b, linestyle='--', alpha=0.5, color='blue')
    plt.text(b+0.1,0.1,str(b),rotation=90, size=8)
plt.ylabel('Kernal density estimate')
plt.xlabel('Average nucleotide identity')
pdf.savefig("ani_histogram.pdf")
plt.close()
