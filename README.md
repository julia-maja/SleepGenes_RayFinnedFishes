# Evolution of visual and non-visual opsins in ray-finned fishes

## I - Opsin gene mining procedure

The mining of opsin genes is done using the script [Opsin_miner.sh](Opsin_miner.sh) 

Every file needed for the script are available in the "Database" folder, except for the uniprot database which can be downloaded here : https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue&facets=reviewed%3Atrue

The script has to be launched like this : 

./[Opsin_miner.sh](Opsin_miner.sh)  Genome_fasta_file.fa [Database/all_opsins_cdhit80.prot](Database/all_opsins_cdhit80.prot) uniprot_sprot.fasta [Database/Scripts_opsins/](Database/Scripts_opsins/) max_intron_size number_of_thread

- Genome_fasta_file.fa = Your genome of interest, in fasta format  
-  max_intron_size = The maximum intron size you want. I used 50000 for all ray-finned fishes genomes.  
- number_of_thread = The number of threads to use. 

Note that the script is intented to be used in a SLURM environment. 

If you work in a non-slum environment, replace those two sbatch commands (line 127 and line 128) : 

`sbatch --job-name=ops_vs_scaff -W -c 2 --qos=6hours --mem=4G --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.exo.rslt ; sleep 10" &`

`sbatch --job-name=ops_vs_scaff -W -c 2 --qos=6hours --mem=4G --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.noexhaustive.exo.rslt ; sleep 10" &`

by = 

`nohup $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.exo.rslt &`

`nohup $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.noexhaustive.exo.rslt &`


If you work under a slurm environment, modify `--qos=6hours` by an existing partition in your environment.



## II - Species tree

The species tree can be generated from BUSCO results using the script [Make_Actino_Species_tree.sh](Make_Actino_Species_tree.sh).

## III - Opsin gene alignments and phylogenies

Gene alignments, phylogenies, gene tree-species tree reconcilations and dN/dS analysis can be reproduced with the script [Opsin_phylogeny_process.sh](Opsin_phylogeny_process.sh). 


## IV - RNA-seq data and mapping

Script to clean RNA-seq data  mapping : [RNA_seq_mapping.sh](RNA_seq_mapping.sh). 

## V - All data analysis

All the analysis presented in the study can be replicated using the different Rscripts provided in this guthub repository. 
All the file needed to run the scripts are available on FigShare, in folders with the same name as the script.

For example, to run the script fish_ospins_ASTRAL.R, you will need the file contained in the folder fish_ospins_ASTRAL in the FigShare repository mentionned above.


[fish_ospins_ASTRAL.R](fish_ospins_ASTRAL.R) : Script to investigate the co-evolution of opsin subfamilies,  the co-evolution between visual opsins and non-visual opsins, and associations between ecological factors and the opsin repertoire

[Rates_Selection_ASTRAL.R](Rates_Selection_ASTRAL.R) : Script to analyse the birth and death rates of opsins, as well as their dN/dS + the association between dN/dS and ecological factors.


[RNAseq_opsins.R](RNAseq_opsins.R) : Script to perform the RNA-seq data anlysis presented in the study

[Robustness_Assembly.R](Robustness_Assembly.R) : Script to compare the opsin repertoire found in different genome assemblies of the same species

[opsins_gene_tree.R](opsins_gene_tree.R) : Script to draw the opsin maximum likelihood tree and the clade tree

[Pinopsin.R](Pinopsin.R) : Script to perform analysis of the pinopsin gene

[opsins_Annot_comparison.R](opsins_Annot_comparison.R) : Script to compare the results of the opsin mining pipeline used in this study and previously reported number of opsins in various ray-finned fishes. It also compare the results with the NCBI annotations.

[Omega_OGG_analysis.R](Omega_OGG_analysis.R) : Script to perform pGLS between the mean dN/dS per orthogroup and the number of opsins, considering all ray-finned fishes

[OMEGA_OGG_pGLS_Teleost.R](OMEGA_OGG_pGLS_Teleost.R) : Script to perform pGLS between the mean dN/dS per orthogroup and the number of opsins, considering only diploid teleost


[pGLS_Opsins_vs_OGG.R](pGLS_Opsins_vs_OGG.R) : Script to perform pGLS between the number of copy number per orthogroup and the number of opsins, considering all ray-finned fishes

[pGLS_Opsins_vs_OGG_teleost.R](pGLS_Opsins_vs_OGG_teleost.R) : Script to perform pGLS between the number of copy number per orthogroup and the number of opsins, considering only diploid teleost

[OGG_analysis_ASTRAL.R](OGG_analysis_ASTRAL.R) : Script to analayse pGLS results between the gene copy number in each orthogroup and the number of opsins.

[GW_Omega_analysis_ASTRAL.R](GW_Omega_analysis_ASTRAL.R) : Script to analayse pGLS results between the mean dN/dS of each orthogroup and the number of opsins.


!! Important : The results of gprofiler (enrichment analysis) in OGG_analysis_ASTRAL.R and GW_Omega_analysis_ASTRAL.R can be inconsistent with the study because of changes in the databases. To have the exact same results as the study, one can use the 2023-09-14 database, which you can access here: https://biit.cs.ut.ee/gprofiler_archive3/e110_eg57_p18/gost !!








