import sys, os, glob, re, shutil
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sb
import numpy as np
from Bio import SeqIO, Seq


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

def euk_pull_proteins(genome_fasta, genome_gff3, protein_output_file):
    ## pulls CDS sequences out of eukaryotic genome gff3 file provided and outputs
    ## protein sequences to use in orthologue calling
    gff_cds = {}
    with open(genome_gff3, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            if line[2] in ['cds', 'CDS']:
                note = line[8].split(";")
                ## pulls out the parent mRNA ID of the cds 
                parent_ind = [note.index(x) for x in note if 'Parent' in x][0]
                parent = re.sub("Parent=", "", note[parent_ind])

                ## pulls out the id of the cds
                id_ind = [note.index(x) for x in note if 'ID=' in x][0]
                name = re.sub("ID=", "", note[id_ind])

                ## grabs the rest of information needed to pull out the CDS sequence from
                ## the fasta file
                chrom = line[0]
                start = int(line[3])
                stop = int(line[4]) 
                strand = line[6]
                #region = '{}:{}-{}'.format(chrom, start, stop)
                if parent not in gff_cds.keys():
                    gff_cds[parent] = [[chrom, start, stop, strand]]
                else:
                    gff_cds[parent].append([chrom, start, stop, strand])

    ## goes through the input fasta file and pulls out CDS sequences
    out_f = open(protein_output_file, 'w')
    record_dict = SeqIO.index(genome_fasta, "fasta")
    for mrna in gff_cds:
        all_cds = gff_cds[mrna]
        full_cds = ''
        start_flag = 0
        for cds in all_cds:
            if start_flag == 0: ## grab the first start
                start = cds[1]
                start_flag = 1
            stop   = cds[2] ## grab the last stop and strand
            strand = cds[3]
            if strand == '-':
                strand = 0
            else:
                strand = 1
            seq = record_dict[cds[0]][cds[1]-1:cds[2]].seq
            if cds[3] == "-":
                seq = seq.reverse_complement()
            full_cds += seq
        while len(full_cds) % 3 > 0: ## adding N's to partial CDS
            full_cds += 'N'
        protein = full_cds.translate()
        outstring = ">{m_id} # {start} # {stop} # {strand}\n{seq}\n".format(m_id = mrna, start = start, stop = stop, strand = strand, seq = protein)
        out_f.write(outstring)
    out_f.close()
    record_dict.close()
    
def run_pyparanoid_pipeline(pyp_dir, genomes, data_dir, threads, gff3 = False):
    ## Set up pyparanoid input directories
    if not os.path.exists(pyp_dir):
        os.mkdir(pyp_dir)
    if not os.path.exists(pyp_dir+'/gdb'):
        os.mkdir(pyp_dir+'/gdb')
    if not os.path.exists(pyp_dir+'/gdb/pep'):
        os.mkdir(pyp_dir+'/gdb/pep')
    if not gff3:
        print("Calling genes (prodigal)")
        with open(pyp_dir+'/strainlist.txt', 'w') as sfh:
            for fna in genomes:    
                prefix = os.path.splitext(os.path.basename(fna))[0]
                ## Call genes with prodigal
                os.system(' '.join(['prodigal', '-t', prefix+'.tmp.ptrain', '-c', '-i', fna, '>', '/dev/null', '2>&1']))
                peptide_out = pyp_dir+'/gdb/pep/'+prefix+'.pep.fa'
                gff3_out = pyp_dir+'/gdb/pep/'+prefix+'.gff3'
                command = 'prodigal -c -i {fasta} -a {pep} -t {train} -f gff -o {gff3} > /dev/null 2>&1'.format(fasta = fna, pep = peptide_out, train = prefix+'.tmp.ptrain', gff3 = gff3_out) 
                os.system(command)
                os.system(' '.join(['cp', fna, pyp_dir+'/gdb/'+prefix+'.fasta']))
                os.system(' '.join(['rm', prefix+'.tmp.ptrain']))
                sfh.write(prefix+"\n")
                rename_prodigal_peptides(peptide_out)
                
    else:
        print("Parsing genes from gff3")
        with open(pyp_dir+'/strainlist.txt', 'w') as sfh:
            for fna in genomes:
                prefix      = os.path.splitext(os.path.basename(fna))[0]
                print(prefix)
                gff3        = [x for x in os.listdir('genomes_gff3') if prefix in x][0]
                gff3        = os.path.join('genomes_gff3', gff3)
                protein_out = pyp_dir+'/gdb/pep/'+prefix+'.pep.fa'
                if os.path.exists(protein_out):
                    pass
                else:
                    ## pull out proteins from gff3 files
                    euk_pull_proteins(genome_fasta = fna,
                                      genome_gff3 = gff3,
                                      protein_output_file = protein_out)
                    os.system(' '.join(['cp', fna, pyp_dir+'/gdb/'+prefix+'.fasta']))
                    sfh.write(prefix+"\n")

    ## copy small mibig fasta files to the orthofinder peptide directory
    src_mibig = os.path.join(data_dir, 'mibig_2.0_diamond')
    src_files = os.listdir(src_mibig)
    for file_name in src_files:
        full_file_name = os.path.join(src_mibig, file_name)
        if os.path.isfile(full_file_name):
            shutil.copy(full_file_name, pyp_dir+'/gdb/pep/')
    ## run orthofinder (the export may not be necessary for all computers, maybe turn that into an exception?
    command = 'orthofinder -t {cpu}  -os -M msa -f {pep}'.format(cpu = str(threads), pep = pyp_dir+'/gdb/pep/')
    print('Running:', command)
    os.system(command)



def parse_orthogroups_tsv(pyp_dir, focal_strain = False):
    ortho_finder_dir = pyp_dir+'/gdb/pep/'+'OrthoFinder'
    out_path = os.path.join(sorted(Path(ortho_finder_dir).iterdir(), key=os.path.getmtime)[0], "Orthogroups")
    all_groups = os.path.join(out_path, 'Orthogroups.tsv')
    ## have to read in with pandas because the whitespace is weird
    df = pd.read_csv(all_groups, sep="\t")
    df.dropna(thresh = 1)
    df.columns = [re.sub('.pep', '', x) for x in df.columns]
    if not focal_strain:
        df_dict    = df.set_index('Orthogroup').T.to_dict('list') ## make it a dictionary
    else:
        focal_strain = os.path.splitext(focal_strain)[0]
        df_sub     = df[["Orthogroup", focal_strain]].dropna() ## only keep values from the focal genome
        df_dict    = df_sub.set_index('Orthogroup').T.to_dict('list') ## make it a dictionary
    return(df_dict)

def parse_pyparanoid_results(pyp_dir, focal_strain, ortholog_tall, ortholog_list, mibig_groups):

    focal_strain = os.path.splitext(focal_strain)[0]
    ortho_finder_dir = pyp_dir+'/gdb/pep/'+'OrthoFinder'
    out_path = os.path.join(sorted(Path(ortho_finder_dir).iterdir(), key=os.path.getmtime)[0], "Orthogroups")

    ## Define input files
    counts     = os.path.join(out_path, 'Orthogroups.GeneCount.tsv')

    ## Define some output files
    mibig_groups = open('mibig_groups.tsv', 'w')
    ortholog_tall = open('ortho.tsv', 'w')
    with open(counts, 'r') as fh:
        line1 = 0
        header = 0
        mb_header = 0
        for line in fh:
            line = line.strip().split()
            if line1 == 0:
                keys = [re.sub('.pep', '', x) for x in line]
                orthogroup_pos = keys.index('Orthogroup')
                keys.pop(orthogroup_pos)
                total_pos      = keys.index('Total')
                keys.pop(total_pos)
                mibig_cols = [keys.index(x) for x in keys if 'sequences' in x] ## pull out the mibig db columns
                keys = keys[:mibig_cols[0]]
                line1 = 1
            else:
                orthogroup   = line.pop(orthogroup_pos)
                total_counts = line.pop(total_pos)
                mibig_count = sum([int(line[x]) for x in mibig_cols])
                line = line[:mibig_cols[0]]
                line_dict    = dict(zip(keys, line))
                line_dict['mibig'] = mibig_count
                ## remove lines that only have mibig hits
                if mibig_count == total_counts:
                    continue
                else:
                    for genome in line_dict:
                        if header == 0:
                            header = "ortho_group\tstrain_id\tn_ortho\n"
                            ortholog_tall.write(header)
                        if genome == 'mibig':
                            continue
                        outstring = "{}\t{}\t{}\n".format(orthogroup, genome, line_dict[genome])
                        ortholog_tall.write(outstring)
                    if mibig_count != '0':
                        if mb_header == 0:
                            mb_header = 'ortho_group\tn_ortho\n'
                            mibig_groups.write(mb_header)
                        outstring = "{}\t{}\n".format(orthogroup, mibig_count)
                        mibig_groups.write(outstring)

    mibig_groups.close()
    ortholog_tall.close()


    ## Define the last output file
    ortholog_list = open('ortho_group_list.tsv', 'w')
    
    df_dict = parse_orthogroups_tsv(pyp_dir = pyp_dir, focal_strain = focal_strain)
    outstring = "{}\t{}\t{}\n"
    ortholog_list.write(outstring.format('strain_id', 'Gene', 'ortho_group'))
    for orthogroup in df_dict:
        gene = df_dict[orthogroup][0]
        if "," in gene:
            genes = re.sub(' ', '', gene).split(",")
            for gene in genes:
                ortholog_list.write(outstring.format(focal_genome, gene, orthogroup))
        else:
            ortholog_list.write(outstring.format(focal_genome, gene, orthogroup))

    ortholog_list.close()





def get_gff_dict(gff3, chr_dict = False):
    if gff3:
        gff_dir = 'genomes_gff3'
        
    else:
        gff_dir = 'pyp/gdb/pep/'
    gff_list = [x for x in os.listdir(gff_dir) if 'gff' in x]
    gff_dict = {}
    for x in gff_list:
        gff_in = os.path.join(gff_dir, x)
        with open(gff_in, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split()
                if line[2] in ['cds', 'CDS']:
                    note = line[8].split(";")
                    genome_id = os.path.splitext(x)[0]
                    chrom = line[0]
                    start = int(line[3])
                    stop = int(line[4]) 
                    strand = line[6]
                    if gff3:
                        ## pulls out the parent mRNA ID of the cds
                        id_ind = [note.index(x) for x in note if 'Parent' in x][0]
                        #name = re.sub("Parent=", "", note[id_ind])
                        gene_id = re.sub("Parent=", "", note[id_ind])
                    else:
                        ## need to make unique for gff3 from prodigal
                        ## pulls out the id of the cds
                        id_ind = [note.index(x) for x in note if 'ID=' in x][0]
                        name = re.sub("ID=", "", note[id_ind])
                        gene_id = genome_id + "-" + chrom + "-" + name
                        
                    ## grabs the rest of information needed to pull out the CDS sequence from
                    ## the fasta file
                    if not chr_dict:
                        if gene_id not in gff_dict.keys():
                            gff_dict[gene_id] = [genome_id, chrom, start, stop, strand]
                        else:
                            gff_dict[gene_id].append([genome_id, chrom, start, stop, strand])
                    else:
                        if genome_id not in gff_dict:
                            gff_dict[genome_id] = {}
                            if chrom not in gff_dict[genome_id]:
                                gff_dict[genome_id][chrom] = [start, stop, strand, gene_id]
                            else:
                                gff_dict[genome_id][genome_id][chrom].append([start, stop, strand, gene_id])
                        else:
                            if chrom not in gff_dict[genome_id]:
                                gff_dict[genome_id][chrom] = [start, stop, strand, gene_id]
                            else:
                                gff_dict[genome_id][chrom].append([start, stop, strand, gene_id])
    return(gff_dict)


def make_gene_map(gff_cds, gene_map):
    ## pulls CDS sequences out of eukaryotic genome gff3 file provided and outputs
    ## protein sequences to use in orthologue calling
    gene_orthos = {}
    orthogroups = parse_orthogroups_tsv(pyp_dir = pyp_dir)
    for i in orthogroups:
        for x in orthogroups[i]:
            try:
                all_genes = [y.strip() for y in x.split(",")]
                for gene in all_genes:
                    if gene in gene_orthos:
                        gene_orthos[gene].append(i)
                    else:
                        gene_orthos[gene] = [i]
            except:
                pass
    ## this removes genes that don't have orthogologues
    with open(gene_map, 'w') as gfh:
        gfh.write("\t".join(['strain_id', 'gene_name', 'contig', 'ortho_group'])+"\n")
        for i in gff_cds:
            if i in gene_orthos.keys():
                if len(gene_orthos[i]) > 1:
                    orthogroup = ','.join(gene_orthos[i])
                else:
                    orthogroup = gene_orthos[i][0]
                gfh.write("\t".join([gff_cds[i][0], gff_cds[i][1], i,  orthogroup]) + "\n")
            else:
                continue



def rename_prodigal_peptides(fasta):
    fasta_temp_name = fasta + '.temp'
    fasta_t = open(fasta_temp_name, 'a')
    genome_prefix   = re.sub(".pep$", "", get_prefix(fasta))
    with open(fasta, 'r') as fh:
        for line in fh:
            if line.startswith(">"):
                line_t = re.sub(">", "", line)
                line_elements = line_t.split("#")
                line_elements = [x.strip() for x in line_elements]
                gff_note      = line_elements[4].split(';')
                gene_id       = re.sub('ID=', '', gff_note[0])
                chrom = re.sub("_[0-9]+.*$", "", line_elements[0]) ## remove prodigal chr extension
                new_id = '>' + genome_prefix + "-" + chrom + "-" + gene_id + "\n"
                fasta_t.write(new_id)
            else:
                fasta_t.write(line)
    os.remove(fasta)
    os.rename(fasta_temp_name, fasta)
                
                




def parse_antismash(antismash_dir, antismash_loci, gene_in_bgc, gff_dict):
    ## Read in antiSMASH BGC genome coordinates
    as5 = {}
    with open(antismash_loci, 'w') as afh:
        for gbk in glob.glob(antismash_dir+"/*.gbk"):
            if not re.search( r"region\d+", gbk ): ## We only want the genome-wide gbk, so ignore orthers
                genome = os.path.splitext(os.path.split(gbk)[1])[0]
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
        gfh.write("\t".join(['genome', 'contig', 'gene_id', 'start', 'end', 'strand', 'inclust'])+"\n")
        for genome in gff_dict:
            if genome in as5:
                for contig in gff_dict[genome]:
                    start = int(gff_dict[genome][contig][0])
                    stop = int(gff_dict[genome][contig][1])
                    strand = gff_dict[genome][contig][2]
                    gene_id = gff_dict[genome][contig][3]
                    inclust = 'False'
                    if contig in as5[genome]:
                        for as_start in as5[genome][contig]:
                            if start >= as_start and start <= as5[genome][contig][as_start]:
                                inclust = 'True'
                            elif stop >= as_start and stop <= as5[genome][contig][as_start]:
                                inclust = 'True'
                    gfh.write("\t".join([genome, contig, gene_id, str(start), str(stop), strand, inclust])+"\n")
            else:
                continue

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




if __name__ == '__main__':
    if len(sys.argv[1:]) < 5:
        sys.stderr.write('''Usage: python %s <focalgenome.fna> <genome_folder> <antismash_results_folder> <ko_list> <prokaryotes.hal> -gff3
        <focalgenome.fna>           nucleotide fasta file of the genome of interest
        <genome_folder>             folder of all genomes (yes, including the focal_genome) to be used as comparators
        <antismash_results_folder>  folder of antiSMASH 5 (or greater) genbank results.
        <ko_list>                   downloaded from ftp://ftp.genome.jp/pub/db/kofam/
        <prokaryotes.hal>           downloaded from ftp://ftp.genome.jp/pub/db/kofam/ with the corresponding HMMs
        <threads>                   the number of threads for parallel computing
        <-gff3>                     flag to tell evomining to look for user supplied annotated genes in directory "genomes_gff3" \n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)
    
    ## Parse command line paramaters
    focal_genome, genome_dir, antismash_dir, ko_list, hal_db = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

    ## set gff3 flag
    if "-gff3" in sys.argv:
        gff3 = True
    else:
        gff3 = False
        
    ## Set internal parameters
    threads = sys.argv[6]
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
        


    
    ## Parse ortholog results
    if os.path.exists(ortholog_tall) and os.path.exists(ortholog_list) and os.path.exists(mibig_groups) and os.path.exists(pyp_dir+'/gdb/pep/OrthoFinder'):
        print('Orthologs already calculated, see results files:', ortholog_tall, ortholog_list, mibig_groups)
    else:
        run_pyparanoid_pipeline(pyp_dir, genomes, data_dir, threads, gff3 = gff3)
        print("Parsing ortholog results...")
        parse_pyparanoid_results(pyp_dir, focal_strain, ortholog_tall, ortholog_list, mibig_groups)

    if not os.path.exists(gene_map):
        ## Map genes
        print("Creating gene map...")
        gff_dict = get_gff_dict(gff3)
        make_gene_map(gene_map = gene_map, gff_cds = gff_dict)

    gff_dict = get_gff_dict(gff3, chr_dict = True)
    ## Parse antiSMASH results
    print("Parsing antiSMASH results...")
    parse_antismash(antismash_dir, antismash_loci, gene_in_bgc, gff_dict)

    ## Functional annotation
    if os.path.exists(kofam_parsed):
        print("Functional annotation already completed and written to "+kofam_parsed)
    else:
        print("Performing functional annotation (kofamscan)...")
        kofam_annotation(ko_list, hal_db, focal_strain, pyp_dir, kofam_out, kofam_parsed)

    ## Identify expansions
    print("Identifying putative metabolic expansions...")
    identify_expansions(focal_strain, data_dir, ani, kofam_parsed, ingroups_file, ortholog_tall, ortholog_list, mibig_groups)


