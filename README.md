# EvoMining3

<center><img src="/data/logos.png" alt="logos"
        title="EvoMining3 is a joint collaboration between the University of Wisconsin-Madison and Cinvestav-LangeBio." width="300" height="200" /></center>

EvoMining3 is a joint collaboration between the University of Wisconsin-Madison and Cinvestav-LangeBio.

## Brief description

EvoMining is an evolution based genome mining tool that performs computational and phylogenomic analyses to identify expansions on central metabolic enzyme families and their recruitment by biosynthetic pathways.

## Set up environment

Dependencies are listed in `conda_environment.yml`. It is highly suggested for users to create their own conda environment using this file, e.g.:

`conda env create --file conda_environment.yml --name evomining3`

This creates a new environment called `evomining3` with all dependencies installed. This environment can now be accessed by:

`conda activate evomining3`

## Running EvoMining3 to predict metabolic expansion events
* `python evomining3_identify_expansions.py <focal_genome.fasta> <genomes_directory> <antismash_gbk_directory> <ko_list> <prokaryote.hal>`
  * `focal_genome.fasta`: nucleotide fasta file of the genome of interest
  * `genomes_directory`: folder of all genomes (yes, including the focal_genome) to be used as comparators
  * `antismash_gbk_directory`: folder of antiSMASH 5 (or greater) genbank results. These are used to define BGC boundaries only, so for the purposes of EvoMining, running antiSMASH with the `--minimal` is sufficient (and will save you a lot of compute time). Of course, the user may have other analyses in mind, but any antiSMASH output, `--minimal` or not is fine. It is assumed that users will generate antiSMASH5 results before EvoMining for each genome in `genomes_directory`. A helper script is included at `src/as500_pipe.pl` that will run antiSMASH5 from a docker on each positional argument. After antiSMASH is run (by whatever means the user finds most convienient) copy all resulting genbank files to a new directory `antismash_gbk_directory`. While EvoMining will only look at the whole genome genbanks (and not the ones for individual BGCs), it is fine to copy those into `antismash_gbk_directory` as EvoMining will ignore any "region" genbank.
  * `ko_list` and `prokaryote.hal` are part of kofamscan and must be downloaded from `ftp://ftp.genome.jp/pub/tools/kofamscan/` with the corresponding HMMs

## Using EvoMining3 to generate a gene neighborhood network to mine for bottlenecks and connected expansions (under development)
* `python src/evomining3_calculate_adjacency_index.py`
* `python src/generate_networks.py`

## Docker developmental stage  
`docker build -t evomarc .`  
`cd evomining3  `  
inputs directory contains the genomes in fasta format  
`docker run -i -t -v $(pwd)/data:/home -v $(pwd)/inputs:/home/fasta evomarc /bin/bash`  
inside de docker we run succesfully ani_fj.pl  
`# ani_fj.pl bacillus_stats_cut.tsv`  
produces an output temporary at home inside docker (data outside docker)  
``
