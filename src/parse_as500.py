#!/bin/env python

import glob, os, re
from Bio import SeqIO


for gbk in glob.glob("as500/*/*/*.gbk"):
    if not re.search( r"region\d+", gbk ):
        genome = os.path.split(gbk)[0].split('/')[-1]
        for seq in SeqIO.parse(gbk, "genbank"):
            for feat in seq.features:
                if feat.type == 'region':
                    print("\t".join([genome, seq.id, str(feat.location.start), str(feat.location.end)]))
