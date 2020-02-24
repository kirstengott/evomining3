# EvoMining3

## Calculate ANI
* `python src/calculate_ani.py <focal_genome.fasta> <genomes_directory>`

## Define subset based on ANI and alignment len
* `python src/identify_ani_groups.py <ANI_file.tsv> <focal_genome.fasta>`
* `python src/split_ani_groups.py ANI.tsv <ANI_file.tsv> <ANI_threshold_1> <ANI_threshold_2>....<ANI_threshold_n>`

## Define ortho groups
* `python src/pyparanoid_pipeline.py <genome_dir>`
* `python create_ortho_tall.py <focal_genome.fasta>`

## antiSMASH
* `perl src/as500_pipe.pl <genome1.fasta> <genome_2.fasta>....<genome_n.fasta>`
  * Note: you can run any antiSMASH version yourself, this is just provided as a helper script
* `python src/parse_antismash.py <pyparanoid_pep_folder> <antismash_gbk_1> <antismash_gbk_2>....<antismash_gbk_n>`

## KOFAM functional annotation
* `python src/kofam_annotation.py <ko_list> <prokaryote.hal> <focal_genome_pep.faa>`
  * Note: ko_list and prokaryote.hal are part of kofamscan and must be downloaded from `ftp://ftp.genome.jp/pub/tools/kofamscan/` with the corresponding HMMs

## Go in again and identify expansions
* src/fj_met_expansions.R

## Adjacency index calculations and network generation
* `perl src/genenamemap.pl > genenamemap.tsv`
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


