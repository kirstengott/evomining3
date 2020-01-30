#!/bin/env python

from __future__ import division
import subprocess, sys, os
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import pygraphviz
import matplotlib.pyplot as plt

def get_ani(ani_file):
    ani = {}
    with open(ani_file, 'r') as anifh:
        next(anifh)
        for ln in anifh.read().splitlines():
            a = ln.split("\t")
            ani[a[0]] = {}
            ani[a[0]] = float(a[1])
    return ani
    
def get_network(group, flank, thresh, ani, ai_file, inclust, soi):
    g = nx.Graph()
    seen_nodes = {}
    ai = subprocess.check_output(" ".join(['grep', group, ai_file]), shell=True).decode("utf-8").split("\n")
    for ln in ai:
        if ln == '':
            continue
        g1, gene1, grp, g2, gene2, f, ai, edge = ln.split("\t")
        if group == grp and flank == int(f):
            if gene1 not in seen_nodes:
                g.add_node(gene1,
                           genome = g1,
                           ani = ani[g1],
                           inclust = inclust[gene1])
                seen_nodes[gene1] = 1
            if gene2 not in seen_nodes:
                g.add_node(gene2,
                           genome = g2,
                           ani = ani[g2],
                           inclust = inclust[gene2])
                seen_nodes[gene2] = 1
            if gene1 != gene2 and float(ai) >= thresh:
                g.add_edge(gene1, gene2,
                           weight = float(ai))

    u = g.copy()
    for comp in nx.connected_components(u):
        has_focal = False
        for node in comp:
            if g.nodes[node]['genome'] == soi:
                has_focal = True
        if not has_focal:
            for node in comp:
                g.remove_node(node)

    return g, u

def plot_network(g, soi, out_file):
    pos = graphviz_layout(g, prog="neato")
    as5_foc, not_as5_foc, as5_rest, not_as5_rest = g.copy(), g.copy(), g.copy(), g.copy()
    for node in g.nodes:
        if g.nodes[node]['genome'] != soi:
            if g.nodes[node]['inclust'] != 'True':
                ## Not focal strain, not in clust
                as5_foc.remove_node(node)
                not_as5_foc.remove_node(node)
                as5_rest.remove_node(node)
            else:
                ## Not focal strain, in clust
                as5_foc.remove_node(node)
                not_as5_foc.remove_node(node)
                not_as5_rest.remove_node(node)
        else:
            if g.nodes[node]['inclust'] != 'True':
                ## Focal strain, not in clust
                as5_rest.remove_node(node)
                not_as5_rest.remove_node(node)
                as5_foc.remove_node(node)
            else:
                ## Focal strain, in clust
                as5_rest.remove_node(node)
                not_as5_rest.remove_node(node)
                not_as5_foc.remove_node(node)
    #all_colors = [ g.nodes[n]['ani'] for n in g]
    au_colors = [ as5_foc.nodes[n]['ani'] for n in as5_foc]
    ar_colors = [ as5_rest.nodes[n]['ani'] for n in as5_rest]
    nu_colors = [ not_as5_foc.nodes[n]['ani'] for n in not_as5_foc]
    nr_colors = [ not_as5_rest.nodes[n]['ani'] for n in not_as5_rest]
    antismash_shape = 'd' # diamond
    other_shape = 'o' # circle
    foc_color = '#d91600' # red outline
    other_color = '#000000' # black outline
    ec = nx.draw_networkx_edges(g, pos,
                                alpha = 0.5,
                                width = 0.1,
                                edge_color = '#555555' # gray
    )
    nr = nx.draw_networkx_nodes(not_as5_rest, pos,
                                node_color = nr_colors,
                                node_shape = other_shape,
                                edgecolors = other_color,
                                linewidths = 0.2,
                                vmin = 85,
                                vmax = 100,
                                node_size = 2,
                                alpha = 0.8,
                                with_labels = False
    )
    ar = nx.draw_networkx_nodes(as5_rest, pos,
                                node_color = ar_colors,
                                node_shape = antismash_shape,
                                edgecolors = other_color,
                                linewidths = 0.2,
                                vmin = 85,
                                vmax = 100,
                                node_size = 2,
                                alpha = 0.8,
                                with_labels = False
    )
    au = nx.draw_networkx_nodes(as5_foc, pos,
                                node_color = au_colors,
                                node_shape = antismash_shape,
                                edgecolors = foc_color,
                                linewidths = 0.2,
                                vmin = 85,
                                vmax = 100,
                                node_size = 2,
                                alpha = 0.8,
                                with_labels = False
    )
    nu = nx.draw_networkx_nodes(not_as5_foc, pos,
                                node_color = nu_colors,
                                node_shape = other_shape,
                                edgecolors = foc_color,
                                linewidths = 0.2,
                                vmin = 85,
                                vmax = 100,
                                node_size = 2,
                                alpha = 0.8,
                                with_labels = False
    )
    try: ## this is to catch plots with no colorscale
        plt.colorbar(nr)
    except:
        print("Error occured in colorbar...no ANI diversity?")
    plt.savefig(out_file, dpi=1600)
    plt.clf()

def describe_components(g, soi):
    c = 1
    comp_dict = {}
    for comp in nx.connected_components(g):
        tot, n, foc_n, as_n = 0, 0, 0, 0
        for node in comp:
            tot += g.nodes[node]['ani']
            n += 1
            if g.nodes[node]['genome'] == soi:
                foc_n += 1
            if g.nodes[node]['inclust'] == 'True':
                as_n += 1
        avg = round(tot/n, 3)
        name = 'component'
        if c < 10:
            name = name+'000'+str(c)
        elif c < 100:
            name = name+'00'+str(c)
        elif c < 1000:
            name = name+'0'+str(c)
        comp_dict[name] = {}
        comp_dict[name]['n'] = n
        comp_dict[name]['foc_n'] = foc_n
        comp_dict[name]['avg'] = avg
        comp_dict[name]['as_n'] = as_n
        c += 1
    return comp_dict


ai_file = './ai_all.tsv'
ani_file = './in_groups.tsv'
ani = get_ani(ani_file)
focal_strain = '89329'
n_to_process = 500

if not os.path.isdir('./networks'):
    os.mkdir('./networks')

grp2chk = {}
with open('./med_expansions.tsv') as efh:
    for ln in efh.read().splitlines():
        ln = ln.replace("\"", "")
        ln = ln.split("\t")
        if ln[0] == 'ortho_group':
            continue
        else:
            grp2chk[ln[0]] = float(ln[3]) ## value = median 100 expansion

inclust = {}
with open('./gene_inclust.tsv', 'r') as ifh:
    for ln in ifh.read().splitlines():
        genome, contig, g_id, start, end, ic = ln.split("\t")
        inclust[g_id] = ic

with open('net_sum.tsv', 'w') as o:
    o.write("\t".join(['ortho_group', 'flank', 'threshold', 'component', 'size', 'n_genes_in_focal_strain', 'n_genes_in_antismash_clusts', 'mean_ani'])+"\n")
    flanks = [2, 5]
    cutoffs = [0.00, 0.25, 0.50, 0.75]
    process = True
    n = 0
    for grp in sorted(grp2chk, key=grp2chk.__getitem__, reverse=True):
        if n == n_to_process:
            break
        if process:
            n += 1
            for f in flanks:
                for thresh in cutoffs:
                    pref = 'networks/'+'_'.join([grp, str(f), str(int(100*thresh))])
                    print("\t".join([grp, str(f), str(thresh), str(grp2chk[grp]), str(n)+' of '+str(n_to_process)]))
                    graph_filt, graph_full = get_network(grp, f, thresh, ani, ai_file, inclust, focal_strain)
                    if not os.path.isfile(pref+"_full.png"):
                        plot_network(graph_full, focal_strain, pref+"_full.png")
                    if not os.path.isfile(pref+"_filt.png"):
                        plot_network(graph_filt, focal_strain, pref+"_filt.png")
                    comp = describe_components(graph_filt, focal_strain)
                    for c in sorted(comp):
                        o.write("\t".join([grp, str(f), str(thresh), c, str(comp[c]['n']), str(comp[c]['foc_n']), str(comp[c]['as_n']), str(comp[c]['avg'])])+"\n")
