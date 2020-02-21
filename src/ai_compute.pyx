#!/bin/env python

from __future__ import division

cpdef get_all_group_members(str group, dict genome_dict):
    g = {}
    for contig in genome_dict:
        for gene_contig_num in genome_dict[contig]:
            if genome_dict[contig][gene_contig_num]['ortho_group'] == group:
                if contig in g.keys():
                    g[contig].append(gene_contig_num)
                else:
                    g[contig] = [gene_contig_num]
    return g

cpdef get_flank_set(dict genome_dict, str contig, int contig_gene_num, int flank_sz):
    lo = contig_gene_num - flank_sz
    hi = contig_gene_num + flank_sz
    pairs = set()
    for i in range(lo, hi+1):
        if i in genome_dict[contig] and i+1 in genome_dict[contig]:
            pairs.add(tuple(sorted([genome_dict[contig][i]['ortho_group'], genome_dict[contig][i+1]['ortho_group']])))
        else: ## flanks extend past the contig
            return set()
    return pairs

cpdef compute_ai(dict gmap, str ref, str q, str ortho_group, int flank):
    ref_dict, q_dict = gmap[ref], gmap[q]
    r_group_members = get_all_group_members(ortho_group, ref_dict)
    q_group_members = get_all_group_members(ortho_group, q_dict)
    for goi_contig in r_group_members:
        for goi_contig_num in r_group_members[goi_contig]:
            ref_set = get_flank_set(ref_dict, goi_contig, goi_contig_num, flank)
            for q_contig in q_group_members:
                for q_contig_num in q_group_members[q_contig]:
                    if ref == q and goi_contig == q_contig and goi_contig_num == q_contig_num:
                        continue
                    q_set = get_flank_set(q_dict, q_contig, q_contig_num, flank)
                    AI = 0.0
                    edge_flag = 'Edge'
                    if len(ref_set) > 0 and len(q_set) > 0:
                        ## intersection / union
                        AI = len(ref_set & q_set) / len(ref_set | q_set)
                        edge_flag = 'NotEdge'
                    if AI > 0:
                        print("\t".join([
                            ref,
                            ref_dict[goi_contig][goi_contig_num]['gene_name'],
                            ortho_group,
                            q,
                            q_dict[q_contig][q_contig_num]['gene_name'],
                            str(flank),
                            str(round(AI, 3)),
                            edge_flag
                        ]))

cpdef ai_wrapper():    
    gmap = {}
    grps = {}
    with open('./genenamemap.tsv') as gfh:
        for ln in gfh.read().splitlines():
            g, gene, c, cgnum, ortho_group = ln.split("\t")
            if g == 'genome': ## header
                continue
            cgnum = int(cgnum)
            if g in gmap.keys():
                if c in gmap[g].keys():
                    gmap[g][c][cgnum] = {}
                    gmap[g][c][cgnum]['gene_name'] = gene
                    gmap[g][c][cgnum]['ortho_group'] = ortho_group
                else:
                    gmap[g][c] = {}
                    gmap[g][c][cgnum] = {}
                    gmap[g][c][cgnum]['gene_name'] = gene
                    gmap[g][c][cgnum]['ortho_group'] = ortho_group
            else:
                gmap[g] = {}
                gmap[g][c] = {}
                gmap[g][c][cgnum] = {}
                gmap[g][c][cgnum]['gene_name'] = gene
                gmap[g][c][cgnum]['ortho_group'] = ortho_group
            if ortho_group in grps:
                if g in grps[ortho_group]:
                    grps[ortho_group][g].append(gene)
                else:
                    grps[ortho_group][g] = [gene]
            else:
                grps[ortho_group] = {}
                grps[ortho_group][g] = [gene]
            
    grp2chk = {}
    with open('./med_expansions.tsv') as efh:
        for ln in efh.read().splitlines():
            ln = ln.replace("\"", "")
            ln = ln.split("\t")
            if ln[0] == 'ortho_group':
                continue
            else:
                grp2chk[ln[0]] = 1

    flanks2compute = [2, 5]
    for grp in sorted(grp2chk):
        seen = {}
        for g1 in sorted(grps[grp]):
            for g2 in sorted(grps[grp]):
                sf = "-".join([g1,g2])
                sr = "-".join([g2,g1])
                if sf in seen:
                    continue
                elif sr in seen:
                    continue
                else:
                    for f in flanks2compute:
                        compute_ai(gmap, g1, g2, grp, f )
                    seen[sf] = 1
                    seen[sr] = 1
