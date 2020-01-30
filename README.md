# EvoMining3

## Calculate ANI
* `perl ani_fj.pl`

## Define subset based on ANI and alignment len
* fj_met_expansions.R
* `perl subset.pl`

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
