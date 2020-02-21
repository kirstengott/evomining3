# EvoMining3

## Calculate ANI
* `python calculate_ani.py <focal_genome.fasta> <genomes_directory>`

## Define subset based on ANI and alignment len
* `python identify_ani_groups.py <ANI_file.tsv> <focal_genome.fasta>`
* `python split_ani_groups.py ANI.tsv <ANI_file.tsv> <ANI_threshold_1> <ANI_threshold_2>....<ANI_threshold_n>`

## Define ortho groups
* `perl pyparanoid_pipe.pl`
* `perl get_group_list.pl > group_list.tsv`
* `perl create_ortho_tall.pl`

## antiSMASH
* `perl as500_pipe.pl pan_flavo/gdb/*.fasta`
* `python parse_as500.py > as500.tsv`
* `python inclust.py > gene_inclust.tsv`

## KOFAM functional annotation
* `perl kofam_pipe.pl pan_flavo/gdb/pep/89329.pep.fa`

## Go in again and identify expansions
* fj_met_expansions.R

## Adjacency index calculations and network generation
* `perl genenamemap.pl > genenamemap.tsv`
* `python cython_setup.py build_ext --inplace`
  * Execute if first time
  * If changes to ai_compute.pyx, execute after delete *.c and *.so files
* `python ai_wrap.py > ai_all.tsv`
* `python generate_networks.py`

## Docker developmental stage  
`docker build -t evomarc .`  
`cd evomining3  `  
inputs directory contains the genomes in fasta format  
`docker run -i -t -v $(pwd)/data:/home -v $(pwd)/inputs:/home/fasta evomarc /bin/bash`  
inside de docker we run succesfully ani_fj.pl  
`# ani_fj.pl bacillus_stats_cut.tsv`  
produces an output temporary at home inside docker (data outside docker)  
``


