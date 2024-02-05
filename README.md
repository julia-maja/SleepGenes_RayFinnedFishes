# Evolution of visual and non-visual opsins in ray-finned fishes

## I - Opsin gene mining procedure

The mining of opsin genes is done using the script [Opsin_miner.sh](Opsin_miner.sh) 

Every file needed for the script are available in the "Database" folder, except for the uniprot database which can be downloaded here : https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue&facets=reviewed%3Atrue

The script has to be launched as this : 

./[Opsin_miner.sh](Opsin_miner.sh)  Genome_fasta_file.fa [Database/all_opsins_cdhit80.prot](Database/all_opsins_cdhit80.prot) uniprot_sprot.fasta [Database/Scripts_opsins/](Database/Scripts_opsins/) max_intron_size number_of_thread

- Genome_fasta_file.fa = Your genome of interest, in fasta format  
-  max_intron_size = The maximum intron size you want. I used 50000 for all ray-finned fishes genomes.  
- number_of_thread = The number of threads to use. 

Note that the script is intented to be used in a SLURM environment. 

## II - Opsin gene alignments and phylogenies
