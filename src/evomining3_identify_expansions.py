import sys, os, glob, re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sb
import numpy as np
from Bio import SeqIO


def is_fasta(genome):
    ## Check if file has a fasta file extension
    ext = os.path.splitext(genome)[-1]
    if ext != '.fasta' and ext != '.fna' and ext != '.fa':
        return False
    else:
        return True

def get_prefix(filename):
    return os.path.splitext(os.path.basename(filename))[0]

def get_genomes(genomes_dir):
    ## Populate genomes list from the genome directory
    genomes = []
    for f in glob.glob(genome_dir+"/*.fasta"):
        genomes.append(f)
    for f in glob.glob(genome_dir+"/*.fna"):
        genomes.append(f)
    for f in glob.glob(genome_dir+"/*.fa"):
        genomes.append(f)
    return genomes

def calculate_ani(focal_genome, focal_strain, genomes, ani_file, threads):
    ## Open output file handle and write header
    with open(ani_file, 'w') as afh:
        afh.write("\t".join(['strain_id', 'focal_genome', 'ANI', 'Mappings', 'TotalFrags', 'QueryAligned'])+"\n")
        ## Loop through each query genome
        for query in genomes:
            ## Perform ANI calculation
            os.system(' '.join(['fastANI', '-q', query, '-r', focal_genome, '--visualize', '-t', str(threads), '-o', 'ani_tmp.txt', '>', '/dev/null', '2>&1']))
            ## Parse results to get the length of ANI alignment
            alilen = 0
            with open('ani_tmp.txt.visual', 'r') as vfh:
                for ln in vfh.read().splitlines():
                    col = ln.split("\t")
                    alilen += abs(float(col[6]) - float(col[7]))+1;
            ## Parse results to get the ANI summary information
            ani, maps, tf = 0, 0, 0
            with open('ani_tmp.txt', 'r') as tfh:
                for ln in tfh.read().splitlines():
                    col = ln.split("\t")
                    ani, maps, tf = col[2], col[3], col[4]
            afh.write("\t".join([get_prefix(query), focal_strain, str(ani), str(maps), str(tf), str(alilen)])+"\n")
            for tmp in ['ani_tmp.txt', 'ani_tmp.txt.visual']:
                os.system(' '.join(['rm', tmp]))
    return genomes

def define_ani_groups(focal_genome, ani_file, ani_pdf, aln_frac):
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
    ani = ani[['strain_id', 'ANI']]
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
    ## Plot the histogram
    hist_plot = sb.distplot(ani['ANI'], bins=40)
    pdf = hist_plot.get_figure()
    for b in breaks:
        plt.axvline(x=b, linestyle='--', alpha=0.5, color='blue')
        plt.text(b+0.1,0.1,str(b),rotation=90, size=8)
    plt.ylabel('Kernal density estimate')
    plt.xlabel('Average nucleotide identity')
    pdf.savefig(ani_pdf)
    plt.close()
    return breaks, ani

def split_ani_groups(focal_strain, ani, breaks, ingroups_file):
    ## Include the focal strain
    fs = pd.DataFrame({ focal_strain: [100] }, columns=['strain_id', 'ANI'])
    ani.append(fs)
    ## Assign groups based on breaks (and add 100 for the focal strain vs everything)
    for threshold in breaks + [100]:
        in_t = ani['ANI'] >= threshold
        ani['in_'+str(threshold)] = in_t
    ## Dump groups to file
    ani.to_csv(ingroups_file, sep="\t")

def run_pyparanoid_pipeline(pyp_dir, genomes, data_dir, threads):
    ## Set up pyparanoid input directories
    os.mkdir(pyp_dir)
    os.mkdir(pyp_dir+'/gdb')
    os.mkdir(pyp_dir+'/gdb/pep')
    with open(pyp_dir+'/strainlist.txt', 'w') as sfh:
        for fna in genomes:    
            prefix = os.path.splitext(os.path.basename(fna))[0]
            ## Call genes with prodigal
            os.system(' '.join(['prodigal', '-t', prefix+'.tmp.ptrain', '-c', '-i', fna, '>', '/dev/null', '2>&1']))
            os.system(' '.join(['prodigal', '-c', '-i', fna, '-a', pyp_dir+'/gdb/pep/'+prefix+'.pep.fa', '-t', prefix+'.tmp.ptrain', '>', '/dev/null', '2>&1']))
            os.system(' '.join(['cp', fna, pyp_dir+'/gdb/'+prefix+'.fasta']))
            os.system(' '.join(['rm', prefix+'.tmp.ptrain']))
            sfh.write(prefix+"\n")
    ## Run pyparanoid
    os.system(' '.join(['BuildGroups.py', '--clean', '--verbose', '--use_MP', '--cpus', str(threads), pyp_dir+'/gdb', pyp_dir+'/strainlist.txt', pyp_dir+'/pyparanoid', '>', '/dev/null', '2>&1']))
    os.system(' '.join(['cp', pyp_dir+'/pyparanoid/strainlist.txt', pyp_dir+'/pyparanoid/prop_strainlist.txt']))
    os.system(' '.join(['cp', pyp_dir+'/pyparanoid/homolog.faa', pyp_dir+'/pyparanoid/prop_homolog.faa']))
    os.system(' '.join(['IdentifyOrthologs.py', '--use_MP', '--threshold', '0.9', pyp_dir+'/pyparanoid', pyp_dir+'/ortho', '>', '/dev/null', '2>&1']))
    ## Propogate to MIBiG
    mibig_dir = data_dir+'/mibig_2.0'
    mibig_list = data_dir+'/mibig_2.0_strainlist.txt'
    os.system(' '.join(['PropagateGroups.py', mibig_dir, mibig_list, pyp_dir+'/pyparanoid', '>', '/dev/null', '2>&1']))

def parse_pyparanoid_results(pyp_dir, focal_strain, ortholog_tall, ortholog_list, mibig_groups):
    ## Define inputs
    homolog_faa = pyp_dir+'/pyparanoid/homolog.faa'
    homolog_mat = pyp_dir+'/pyparanoid/homolog_matrix.txt'
    ## Generate the ortholog group list for the focal genome
    with open(ortholog_list, 'w') as gfh:
        gfh.write("\t".join(['strain_id', 'Gene', 'ortho_group'])+"\n")
        with open(homolog_faa, 'r') as ffh:
            for ln in ffh.read().splitlines():
                if ln.startswith('>'):
                    sid, gene, og = ln[1:].split('|')
                    if sid == focal_strain:
                        gfh.write("\t".join([sid, gene, og])+"\n")
    ## Read in the homolog matrix and generate the tall ortholog table
    header = []
    mibig = {}
    with open(ortholog_tall, 'w') as ofh:
        ofh.write("\t".join(['ortho_group', 'strain_id', 'n_ortho'])+"\n")
        with open(homolog_mat, 'r') as mfh:
            for ln in mfh.read().splitlines():
                if len(header) == 0:
                    header = ln[1:].split("\t")
                else:
                    seen = {}
                    group, counts = ln.split("\t", 1)
                    for i, count in enumerate(counts.split("\t")):
                        if header[i].startswith('BGC'):
                            if float(count) > 0:
                                if group in mibig:
                                    mibig[group] += 1
                                else:
                                    mibig[group] = 1
                        else:
                            if header[i] not in seen:
                                ofh.write("\t".join([group, header[i], str(float(count)/2)])+"\n")
                                seen[header[i]] = 1
    with open(mibig_groups, 'w') as mgfh:
        mgfh.write("\t".join(['ortho_group', 'n_ortho'])+"\n")
        for group in mibig:
            mgfh.write("\t".join([group, str(mibig[group])])+"\n")

def parse_antismash(antismash_dir, antismash_loci, gene_in_bgc):
    ## Read in antiSMASH BGC genome coordinates
    as5 = {}
    with open(antismash_loci, 'w') as afh:
        for gbk in glob.glob(antismash_dir+"/*.gbk"):
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
    with open(gene_in_bgc, 'w') as gfh:
        gfh.write("\t".join(['genome', 'contig', 'gene_id', 'start', 'end', 'inclust'])+"\n")
        for faa in glob.glob("pyp/gdb/pep/*.pep.fa"):
            genome = os.path.split(faa)[-1].split('.')[0]
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
    return inclust

def kofam_annotation(ko_list, hal_db, focal_strain, pyp_dir, kofam_out, kofam_parsed):
    ## Note: this is part of kofamscan and a *hal database must be downloaded from
    ## ftp://ftp.genome.jp/pub/tools/kofamscan/
    ## with the corresponding *hmms
    focal_pep = pyp_dir+'/gdb/pep/'+focal_strain+'.pep.fa'
    ## Functional annotation with kofamscan
    if not os.path.isfile(kofam_out):
        os.system(' '.join(['exec_annotation', '-k', ko_list, '-p', hal_db, '--cpu='+str(threads), '-o', kofam_out, focal_pep]))
    ## Parse the results
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
                    kfh.write("\t".join([focal_strain, gene, ko])+"\n")

def identify_expansions(focal_strain, data_dir, ani, kofam_parsed, ingroups_file, ortholog_tall, ortholog_list, mibig_groups):
    ## Parse group and ortho files
    i = pd.read_csv(ingroups_file, delimiter="\t")
    ot = pd.read_csv(ortholog_tall, delimiter="\t")
    ## Left join with in_groups and subset based on ANI criteria
    ot = pd.merge(ot, i, 'left', on='strain_id')
    ot = ot[ot.strain_id.isin(ani.strain_id)]
    ## Get functional data
    func = pd.read_csv(kofam_parsed, delimiter="\t")
    func = func[['Gene', 'KO-D']]
    gl = pd.read_csv(ortholog_list, delimiter="\t")
    gl = gl[['Gene', 'ortho_group']]
    func = pd.merge(func, gl, 'left', on='Gene')
    kegg_hier = data_dir+'/kegg_hier.tsv'
    kh = pd.read_csv(kegg_hier, delimiter="\t")
    func = pd.merge(func, kh, 'left', on='KO-D')
    ## Check if ortho_groups are in metabolism and/or MIBiG
    met_a = ['Metabolism'] ## Kegg level A matches for metabolism
    vitaa_b = ['Amino acid metabolism', 'Metabolism of other amino acids', 'Metabolism of cofactors and vitamins'] ## Kegg level B matches for amino acid or vitamin metabolism
    mb = pd.read_csv(mibig_groups, delimiter="\t")        
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
    with open(exp_list, 'w') as efh:
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
        pdf = hist_plot.get_figure()
        pdf.axes[0].set_yscale('log')
        pdf.savefig(grp+'.hist.pdf')
        plt.close()
    return i

def make_gene_map(gene_map, pyp_dir):
    homolog_fasta = pyp_dir+'/pyparanoid/homolog.faa'
    with open(gene_map, 'w') as gfh:
        gfh.write("\t".join(['strain_id', 'gene_name', 'contig', 'contig_gene_num', 'ortho_group'])+"\n")
        with open(homolog_fasta, 'r') as hfh:
            for ln in hfh.read().splitlines():
                if ln.startswith('>'):
                    genome, gene_name, ortho_group = ln[1:].split("|")
                    contig = re.sub("_\d+$", '', gene_name)
                    contig_gene_num = re.sub(".+_(\d+)$", '\\1', gene_name)
                    gfh.write("\t".join([genome, gene_name, contig, contig_gene_num, ortho_group])+"\n")



if __name__ == '__main__':
    if len(sys.argv[1:]) < 5:
        sys.stderr.write('''Usage: python %s <focalgenome.fna> <genome_folder> <antismash_results_folder> <ko_list> <prokaryotes.hal>
        <focalgenome.fna>           nucleotide fasta file of the genome of interest
        <genome_folder>             folder of all genomes (yes, including the focal_genome) to be used as comparators
        <antismash_results_folder>  folder of antiSMASH 5 (or greater) genbank results.
        <ko_list>                   downloaded from ftp://ftp.genome.jp/pub/db/kofam/
        <prokaryotes.hal>           downloaded from ftp://ftp.genome.jp/pub/db/kofam/ with the corresponding HMMs\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)
    
    ## Parse command line paramaters
    focal_genome, genome_dir, antismash_dir, ko_list, hal_db = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

    ## Set internal parameters
    threads = 60
    aln_frac = 0.25 ## Fraction of alignment cutoff: (minimum number of bp in alignment) / (bp in focal genome)

    ## Check focal genome is correct file extension
    if not is_fasta(focal_genome):
        raise Exception("focal_genome does not have proper extension (.fasta, .fna, or .fa)\n")

    ## Get focal genome name and genomes list
    focal_strain = get_prefix(focal_genome)
    genomes = get_genomes(genome_dir)
    data_dir = os.path.realpath(__file__).replace('src/' + os.path.basename(__file__), 'data')

    
    ## Define output files and folders
    ani_file = 'ANI.tsv'
    ani_pdf = 'ANI_histogram.pdf'
    ingroups_file = 'in_groups.tsv'
    pyp_dir = 'pyp'
    ortholog_tall = 'ortho.tsv'
    ortholog_list = 'ortho_group_list.tsv'
    mibig_groups = 'mibig_groups.tsv'
    antismash_loci = 'antismash.tsv'
    gene_in_bgc = 'gene_in_bgc.tsv'
    kofam_out = focal_strain+'.kofam.txt'
    kofam_parsed = focal_strain+'.kofam.parsed.tsv'
    exp_list = 'expansions_list.tsv'
    gene_map = 'gene_map.tsv'

    ## Calculate ANI
    if os.path.isfile(ani_file):
        print("ANI already calculated and written to "+ani_file)
    else:
        print("Calculating ANI against "+focal_strain+"...")
        calculate_ani(focal_genome, focal_strain, genomes, ani_file, threads)

    ## Define ANI groups
    print("Calculating ANI frequency distribution...")
    breaks, ani = define_ani_groups(focal_genome, ani_file, ani_pdf, aln_frac)
    print("Breakpoints to define ANI groups are:\t"+str(breaks))

    ## Split based on these group breakpoints
    print("Splitting into groups...")
    split_ani_groups(focal_strain, ani, breaks, ingroups_file)

    ## Define orthologs
    if os.path.exists(pyp_dir):
        print("Orthologs already calculated within directory "+pyp_dir)
    else:
        print("Calling genes (prodigal) and calculating orthologs (pyparanoid)...")
        run_pyparanoid_pipeline(pyp_dir, genomes, data_dir, threads)

    ## Parse ortholog results
    print("Parsing ortholog results...")
    parse_pyparanoid_results(pyp_dir, focal_strain, ortholog_tall, ortholog_list, mibig_groups)

    ## Parse antiSMASH results
    print("Parsing antiSMASH results...")
    inclust = parse_antismash(antismash_dir, antismash_loci, gene_in_bgc)

    ## Functional annotation
    if os.path.exists(kofam_parsed):
        print("Functional annotation already completed and written to "+kofam_parsed)
    else:
        print("Performing functional annotation (kofamscan)...")
        kofam_annotation(ko_list, hal_db, focal_strain, pyp_dir, kofam_out, kofam_parsed)

    ## Identify expansions
    print("Identifying putative metabolic expansions...")
    identify_expansions(focal_strain, data_dir, ani, kofam_parsed, ingroups_file, ortholog_tall, ortholog_list, mibig_groups)

    ## Map genes
    print("Creating gene map...")
    make_gene_map(gene_map, pyp_dir)
