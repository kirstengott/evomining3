# EvoMining3

<center><img src="/data/logos.png" alt="logos"
        title="EvoMining3 is a joint collaboration between the University of Wisconsin-Madison and Cinvestav-LangeBio." width="300" height="200" /></center>

EvoMining3 is a joint collaboration between the University of Wisconsin-Madison and Cinvestav-LangeBio.

## Usage
* `python evomining3.py <focal_genome.fasta> <genomes_directory> <antismash_gbk_directory> <ko_list> <prokaryote.hal>`
  * `focal_genome.fasta`: nucleotide fasta file of the genome of interest
  * `genomes_directory`: folder of all genomes (yes, including the focal_genome) to be used as comparators
  * `antismash_gbk_directory`: folder of antiSMASH 5 (or greater) genbank results. It is assumed that users will generate this before EvoMining. A helper script is included at `src/as500_pipe.pl` that will run antiSMASH5 from a docker on each positional argument. After antiSMASH is run (by whatever means the user finds most convienient) copy all resulting genbank files to a new directory `antismash_gbk_directory`. While EvoMining will only look at the whole genome genbanks (and not the ones for individual BGCs), it is fine to copy those into `antismash_gbk_directory` as EvoMining will ignore any "region" genbank.
  * `ko_list` and `prokaryote.hal` are part of kofamscan and must be downloaded from `ftp://ftp.genome.jp/pub/tools/kofamscan/` with the corresponding HMMs

## Still in development...adjacency index calculations and network generation
* `python src/cython_setup.py build_ext --inplace`
  * Execute if first time
  * If changes to ai_compute.pyx, execute after delete *.c and *.so files
* `python src/ai_wrap.py > ai_all.tsv`
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