#Extract all BUSCO genes present in at-least one species

sed -i 's/\.//g' BUSCO90teleost_BUSCO80nonteleost_species.txt #BUSCO90teleost_BUSCO80nonteleost_species.txt = list of species of the dataset
mkdir Common_busco_genes
for species in `cat BUSCO90teleost_BUSCO80nonteleost_species.txt` ; do
	ls -l BUSCO_$species/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/ | sed 's/.* //g' > Common_busco_genes/$species.singlecopy.list
done


cat Common_busco_genes/*.list | grep ".faa" | sort | uniq > all_busco_in_my_species.txt 
cat Common_busco_genes/*.list | grep ".faa" | sort | uniq -c > BUSCO_nb_observations.txt 

sed -i 's/^ *//g' BUSCO_nb_observations.txt  
sed -i 's/ /,/g' BUSCO_nb_observations.txt  


#Extract the 2000 most represented genes

sort -k1 -n -t "," BUSCO_nb_observations.txt | tail -2000 | cut -f2 -d "," > Most_represented_BUSCO_2000.txt


#Align those genes

mkdir Alignements
for i in `cat Most_represented_BUSCO_2000.txt` ; do 
	for species in `cat BUSCO90teleost_BUSCO80nonteleost_species.txt` ; do
		if grep -q "^$i$" Common_busco_genes/$species.singlecopy.list ; then
			sed "s/>.*/>$species/g" ../ALL_Genomes/BUSCO_$species/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/$i >> Alignements/$i
		fi
	done
done


cd Alignements/

for i in *.faa ; do 
	muscle5.1.linux_intel64 -align $i -output $i.aln
done


#Trim all the alignments
for i in *.aln ; do 
	trimal -in $i -out $i.trimal -automated1
done
sed -i 's/?/-/g' *.trimal


#Infer the best model of each gene (G4) and the maximum likelihood tree

for gene_alignment in *.trimal ; do 
	iqtree -s $gene_alignment --seqtype AA -m TEST --mrate G4 -nt 10 -bb 1000
done


# Collapse every nodes that have a low bootstrap value


for gene_tree in *.treefile ; do 
	./collapse_lowsupport.sh $gene_tree $gene_tree.collapsed
done




#Concatenate all trimed alignments

cd ../

python3 AMAS.py concat -f fasta -d aa -i Alignements/*.trimal --part-format nexus

sed -i 's/?/-/g' concatenated.out
mv concatenated.out AMAS_concatenated_alignment_2000BUSCO.fa
mv partitions.txt AMAS_concatenated_alignment_2000BUSCO.partition.nexus

#Create a file with the calibration points

#nano Actino_calibration.txt : 
#Danio_rerio,Megalops_atlanticus	-250
#Danio_rerio,Polypterus_senegalus	-396
#Denticeps_clupeoides,Gasterosteus_aculeatus_aculeatus	-224



#concatenate all gene trees (branch with low support collapsed)
cat Alignements/*.collapsed > All_genes_tree.concatnwk 

java -jar astral.5.7.8.jar -i All_genes_tree.concatnwk -o Astral_nolength.nwk



#Run those R commands to parse the ASTRAL result
>R 
>library('ape')
>library("phytools")
>mytree <- read.tree("Astral_nolength.nwk")
>unrooted_tree <- unroot(mytree)
>write.tree(unrooted_tree, "Astral_unrooted_nolength.nwk")
>outgroup_MRCA <- findMRCA(unrooted_tree, tips=c("Polypterus_senegalus","Erpetoichthys_calabaricus"), type="node")
>mytree_rooted <- root(unrooted_tree, node=outgroup_MRCA, resolve.root= TRUE)
>write.tree(mytree_rooted, "Astral_rooted_nolength.nwk")


sed -i 's/NaN/0.01/g' Astral_unrooted_nolength.nwk
sed -i 's/NaN/0.01/g' Astral_rooted_nolength.nwk


#Date the tree with the least-square method and root with polypteriformes
iqtree2 -s AMAS_concatenated_alignment_2000BUSCO.fa --date Actino_calibration.txt -te Astral_unrooted_nolength.nwk --date-options "-u 0.1" --date-tip 0 --date-ci 100 -o "Polypterus_senegalus,Erpetoichthys_calabaricus" -nt 30 -m LG+F+G4


## Add node labels to the ASTRAL tree using the R commands below

>R
>library("ape")
>mytree <- read.tree("AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk")
>mytree_nodelabel <- makeNodeLabel(mytree, method = "number", prefix = "Node")
>write.tree(mytree_nodelabel, file="AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel")






#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
### Accesstory scripts ################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################





============================================ collapse_lowsupport.sh 


#!/bin/bash


original_tree=$1
collapsed_tree=$2


Rscript collapse_lowsupport.R $original_tree $collapsed_tree




============================================ collapse_lowsupport.R


library(ape)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
Badnodes <- which(as.numeric(mytree$node.label) <= 10) + length(mytree$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){ Badnodes_indexes <- c(Badnodes_indexes, which(mytree$edge[,2] == node)) }

mytree$edge.length[Badnodes_indexes] <- 0 
tree_multi <- di2multi(mytree) 
write.tree(tree_multi, file = args[2])



