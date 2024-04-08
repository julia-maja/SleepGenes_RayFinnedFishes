# Start with these 3 files:
#Pseudogene_opsins.fa => 2207 sequences
#Incomplete_opsins.fa => 1090 sequences
#Complete_opsins.fa => 17318 sequences


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################# Make a phylogeny with every opsin genes !! ######################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

transeq Complete_opsins.fa Complete_opsins.prot ; sed -i 's/_1$//g' Complete_opsins.prot
makeblastdb -in Complete_opsins.prot -dbtype prot
makeblastdb -in Complete_opsins.fa -dbtype nucl


#File with outgroup sequences from sharks and lampreys => Chondrychties_outgroups.prot

mkdir ALL_opsins_phylogeny

cd ALL_opsins_phylogeny/
cp ../Complete_opsins.prot ./

grep ">Lepisosteus_oculatus" Complete_opsins.prot | sed 's/>//g' | sort > template.id
xargs samtools faidx Complete_opsins.prot < template.id > template.prot
grep ">" Complete_opsins.prot | sed 's/>//g' | sort > Complete_opsins.id
comm -23 Complete_opsins.id template.id > rest_of_sequences


xargs samtools faidx Complete_opsins.prot < rest_of_sequences > rest_of_sequences.prot
cat ../Chondrychties_outgroups.prot >> rest_of_sequences.prot
grep ">" ../Chondrychties_outgroups.prot  | sed 's/>//g'  > sharks_id_seq

/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align template.prot -output template.aln
trimal -in template.aln -out template.aln.trimal -gt 0.5 -cons 60

mafft --add rest_of_sequences.prot --keeplength template.aln.trimal > Complete_opsins_plus_outgroup.mafft.aln

sbatch --qos=1day -c 40 --mem=50G iqtree_all_phylo.sh Complete_opsins_plus_outgroup.mafft.aln


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################# Make phylogeny per subfamily ######################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


########## PINOPSIN ######## 

#For pinopsin, no need for trimal step, as the number of sequences is very low and the alignment looks incredibly nice without trimming
mkdir pinopsin
cp ../Clades/pinopsin ./
xargs samtools faidx ../Complete_opsins.prot < pinopsin > pinopsin.prot
xargs samtools faidx ../Complete_opsins.fa < pinopsin > pinopsin.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007894735.2_Callorhinchus-milii-pinopsin >> pinopsin.prot
samtools faidx ../Chondrychties_outgroups.prot XP_041052689.1_Carcharodon-carcharias-pinopsin >> pinopsin.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align pinopsin.prot -output pinopsin.aln
iqtree -s pinopsin.aln --seqtype AA -nt 8 -bb 1000 
echo "XP_007894735.2_Callorhinchus-milii-pinopsin" > outgroup.id
echo "XP_041052689.1_Carcharodon-carcharias-pinopsin" >> outgroup.id
cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R pinopsin.aln.treefile outgroup.id pinopsin.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R pinopsin.aln.nooutgroup.treefile pinopsin.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R pinopsin.aln.nooutgroup.treefile pinopsin.aln.nooutgroup.collapsed95.treefile
xargs samtools faidx pinopsin.aln < pinopsin > pinopsin.fish.aln
/scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/pal2nal.v14/pal2nal.pl pinopsin.fish.aln pinopsin.fa -output fasta > pinopsin.cds.aln
cp pinopsin.cds.aln pinopsin.fish.cds.trimmed.aln

cp pinopsin.aln.nooutgroup.collapsed90.treefile pinopsin.prot.trimmed.aln.nooutgroup.collapsed90.treefile
cp pinopsin.aln.nooutgroup.collapsed95.treefile pinopsin.prot.trimmed.aln.nooutgroup.collapsed95.treefile
cp pinopsin.aln.nooutgroup.treefile parapinopsin.prot.trimmed.aln.nooutgroup.treefile



cd ../



########## Parapinopsin ######## 

mkdir parapinopsin ; cd parapinopsin
cp ../Clades/parapinopsin ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < parapinopsin > parapinopsin.prot
xargs samtools faidx ../Complete_opsins.fa < parapinopsin > parapinopsin.fa


grep "Polypterus\|calabaricus\|Lepiso\|calva\|Danio_rerio" ../Clades/parapinopsin > random_id



sort random_id > temp ; mv temp random_id
grep ">" parapinopsin.fa | sed 's/>//g' | sort > all.id
comm -23 all.id random_id > rest_of_sequences
xargs samtools faidx parapinopsin.prot < random_id > random.prot
xargs samtools faidx parapinopsin.prot < rest_of_sequences > rest_of_sequences.prot
samtools faidx ../Chondrychties_outgroups.prot V9KWA1_Callorhinchus-milii-parapinopsin >> rest_of_sequences.prot

/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align random.prot -output random.aln

mafft --add rest_of_sequences.prot --keeplength random.aln > parapinopsin.aln

xargs samtools faidx parapinopsin.aln < parapinopsin > parapinopsin.fish.aln

trimal -in parapinopsin.fish.aln -backtrans parapinopsin.fa -automated1 -out parapinopsin.fish.cds.trimmed.aln
cp parapinopsin.aln  parapinopsin.prot.trimmed.aln

sbatch --qos=6hours -c 20 --mem=50G iqtree_phylo.sh parapinopsin.prot.trimmed.aln

echo "V9KWA1_Callorhinchus-milii-parapinopsin" > outgroup.id

xargs samtools faidx parapinopsin.prot.trimmed.aln < parapinopsin > parapinopsin.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R parapinopsin.prot.trimmed.aln.treefile outgroup.id parapinopsin.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R parapinopsin.prot.trimmed.aln.nooutgroup.treefile parapinopsin.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R parapinopsin.prot.trimmed.aln.nooutgroup.treefile parapinopsin.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## Opn4x ######## 
 

mkdir opn4x ; cd opn4x
cp ../Clades/opn4x ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn4x > opn4x.prot
xargs samtools faidx ../Complete_opsins.fa < opn4x > opn4x.fa
samtools faidx ../Chondrychties_outgroups.prot NP_001279400.1_Callorhinchus-milii-opn4x >> opn4x.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn4x.prot -output opn4x.aln

xargs samtools faidx opn4x.aln < opn4x > opn4x.fish.aln

trimal -in opn4x.fish.aln -backtrans opn4x.fa -automated1 -out opn4x.fish.cds.trimmed.aln
trimal -in opn4x.aln -automated1 -out opn4x.prot.trimmed.aln

sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn4x.prot.trimmed.aln

echo "NP_001279400.1_Callorhinchus-milii-opn4x" > outgroup.id

xargs samtools faidx opn4x.prot.trimmed.aln < opn4x > opn4x.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn4x.prot.trimmed.aln.treefile outgroup.id opn4x.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn4x.prot.trimmed.aln.nooutgroup.treefile opn4x.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn4x.prot.trimmed.aln.nooutgroup.treefile opn4x.prot.trimmed.aln.nooutgroup.collapsed95.treefile


########## Opn4m2 ######## 


mkdir opn4m2 ; cd opn4m2
cp ../Clades/opn4m2 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn4m2 > opn4m2.prot
xargs samtools faidx ../Complete_opsins.fa < opn4m2 > opn4m2.fa
samtools faidx ../Chondrychties_outgroups.prot NP_001279357.1_Callorhinchus-milii-opn4m1_3 >> opn4m2.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn4m2.prot -output opn4m2.aln
xargs samtools faidx opn4m2.aln < opn4m2 > opn4m2.fish.aln
trimal -in opn4m2.fish.aln -backtrans opn4m2.fa -automated1 -out opn4m2.fish.cds.trimmed.aln
trimal -in opn4m2.aln -automated1 -out opn4m2.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn4m2.prot.trimmed.aln
echo "NP_001279357.1_Callorhinchus-milii-opn4m1_3" > outgroup.id
xargs samtools faidx opn4m2.prot.trimmed.aln < opn4m2 > opn4m2.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn4m2.prot.trimmed.aln.treefile outgroup.id opn4m2.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn4m2.prot.trimmed.aln.nooutgroup.treefile opn4m2.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn4m2.prot.trimmed.aln.nooutgroup.treefile opn4m2.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## Sws1 ######## 

#>Geotria_australis_AY366495_sws1
#MSGDEEFYLFKNISKVGPWDGPQFHIAPKWAFYLQAAFMGFVFICGTPLNAIVLVVTIKYKKLRQPLNYILVNISAAGLVFCLFSISTVFVASMQGYFFLGPTICALEAFFGSLAGLVTGWSLAFLAAERYIVICKPFGNFRFGSKHALVAVGLTWMLGLSVALPPFFGWSRYIPEGLQCSCGPDWYTVGTKYKSEYYTYFLFVFCFVVPLSIIIFSYGSLLGTLRAVAAQQQESASTQKAEREVSRMVIMMVASFCTCYVPYAALAVYMVTNRDHNIDLRFVTVPAFFSKASCVYNPLIYSFMNKQFRACILETVCGKPITDESETSSSRTEVSSVSTTQMIPG


mkdir sws1 ; cd sws1
cp ../Clades/sws1 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < sws1 > sws1.prot
xargs samtools faidx ../Complete_opsins.fa < sws1 > sws1.fa
cat ../sws1_lamprey.prot >> sws1.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align sws1.prot -output sws1.aln
xargs samtools faidx sws1.aln < sws1 > sws1.fish.aln
trimal -in sws1.fish.aln -backtrans sws1.fa -automated1 -out sws1.fish.cds.trimmed.aln
trimal -in sws1.aln -automated1 -out sws1.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh sws1.prot.trimmed.aln
echo "Geotria_australis_AY366495_sws1" > outgroup.id
xargs samtools faidx sws1.prot.trimmed.aln < sws1 > sws1.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R sws1.prot.trimmed.aln.treefile outgroup.id sws1.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R sws1.prot.trimmed.aln.nooutgroup.treefile sws1.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R sws1.prot.trimmed.aln.nooutgroup.treefile sws1.prot.trimmed.aln.nooutgroup.collapsed95.treefile


########## Parietopsin ######## 

#>Lethenteron_camtschaticum_LC600861_parietopsin
#MQEQSGSTAVAVPGAVTPGPTIFPPAGYSFLAAIMFLDASLSIVNNTLVIVITCRYPSLRSPLNALILSMCVSDLLMSLCGTTIAMVSNFHGSLRHIGHAGCVFQGFSVNYFGCVSLWSLTLLAYERRLVVTHGCSMRGGWPRARRGLAFVWTFCLVWAVAPLLGWSAYGPEGVQTSCSIAWERRSLSNYTYLVSYFLACFVIPVSIIVFSYGNVLCSLHTLNKKIKRVGGHPDPREEMRAAVMVLAMVGAFMACWLPYTVLALCVVLSPGTQIPPLVATLPMYFAKTSPIYNPIIYFFLNRQFRACAVEFVTCGAVKLNEPKDESAAPPLPADTAEPPGPTPPRHNQVSPA

mkdir parietopsin ; cd parietopsin
cp ../Clades/parietopsin ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < parietopsin > parietopsin.prot
xargs samtools faidx ../Complete_opsins.fa < parietopsin > parietopsin.fa
cat ../parietopsin_lamprey.prot >> parietopsin.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align parietopsin.prot -output parietopsin.aln
xargs samtools faidx parietopsin.aln < parietopsin > parietopsin.fish.aln
trimal -in parietopsin.fish.aln -backtrans parietopsin.fa -automated1 -out parietopsin.fish.cds.trimmed.aln
trimal -in parietopsin.aln -automated1 -out parietopsin.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh parietopsin.prot.trimmed.aln
echo "Lethenteron_camtschaticum_LC600861_parietopsin" > outgroup.id
xargs samtools faidx parietopsin.prot.trimmed.aln < parietopsin > parietopsin.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R parietopsin.prot.trimmed.aln.treefile outgroup.id parietopsin.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R parietopsin.prot.trimmed.aln.nooutgroup.treefile parietopsin.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R parietopsin.prot.trimmed.aln.nooutgroup.treefile parietopsin.prot.trimmed.aln.nooutgroup.collapsed95.treefile


########## Exorh ######## 


mkdir exorh ; cd exorh
cp ../Clades/exorh ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < exorh > exorh.prot
xargs samtools faidx ../Complete_opsins.fa < exorh > exorh.fa
samtools faidx ../Chondrychties_outgroups.prot A0A7T8R2L6_Rhincodon-typus-rh1 >> exorh.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align exorh.prot -output exorh.aln
xargs samtools faidx exorh.aln < exorh > exorh.fish.aln
trimal -in exorh.fish.aln -backtrans exorh.fa -automated1 -out exorh.fish.cds.trimmed.aln
trimal -in exorh.aln -automated1 -out exorh.prot.trimmed.aln
sbatch --qos=6hours -c 20 --mem=50G iqtree_phylo.sh exorh.prot.trimmed.aln
echo "A0A7T8R2L6_Rhincodon-typus-rh1" > outgroup.id
xargs samtools faidx exorh.prot.trimmed.aln < exorh > exorh.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R exorh.prot.trimmed.aln.treefile outgroup.id exorh.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R exorh.prot.trimmed.aln.nooutgroup.treefile exorh.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R exorh.prot.trimmed.aln.nooutgroup.treefile exorh.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## Tmt1 ########


mkdir tmt1 ; cd tmt1
cp ../Clades/tmt1 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < tmt1 > tmt1.prot
xargs samtools faidx ../Complete_opsins.fa < tmt1 > tmt1.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007892895.2_Callorhinchus-milii-tmt2 >> tmt1.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align tmt1.prot -output tmt1.aln
xargs samtools faidx tmt1.aln < tmt1 > tmt1.fish.aln
trimal -in tmt1.fish.aln -backtrans tmt1.fa -automated1 -out tmt1.fish.cds.trimmed.aln
trimal -in tmt1.aln -automated1 -out tmt1.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh tmt1.prot.trimmed.aln
echo "XP_007892895.2_Callorhinchus-milii-tmt2" > outgroup.id
xargs samtools faidx tmt1.prot.trimmed.aln < tmt1 > tmt1.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R tmt1.prot.trimmed.aln.treefile outgroup.id tmt1.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R tmt1.prot.trimmed.aln.nooutgroup.treefile tmt1.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R tmt1.prot.trimmed.aln.nooutgroup.treefile tmt1.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## tmt2 ########


mkdir tmt2 ; cd tmt2
cp ../Clades/tmt2 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < tmt2 > tmt2.prot
xargs samtools faidx ../Complete_opsins.fa < tmt2 > tmt2.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007892895.2_Callorhinchus-milii-tmt2 >> tmt2.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align tmt2.prot -output tmt2.aln
xargs samtools faidx tmt2.aln < tmt2 > tmt2.fish.aln
trimal -in tmt2.fish.aln -backtrans tmt2.fa -automated1 -out tmt2.fish.cds.trimmed.aln
trimal -in tmt2.aln -automated1 -out tmt2.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh tmt2.prot.trimmed.aln
echo "XP_007892895.2_Callorhinchus-milii-tmt2" > outgroup.id
xargs samtools faidx tmt2.prot.trimmed.aln < tmt2 > tmt2.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R tmt2.prot.trimmed.aln.treefile outgroup.id tmt2.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R tmt2.prot.trimmed.aln.nooutgroup.treefile tmt2.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R tmt2.prot.trimmed.aln.nooutgroup.treefile tmt2.prot.trimmed.aln.nooutgroup.collapsed95.treefile


########## tmt3 ########


mkdir tmt3 ; cd tmt3
cp ../Clades/tmt3 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < tmt3 > tmt3.prot
xargs samtools faidx ../Complete_opsins.fa < tmt3 > tmt3.fa
samtools faidx ../Chondrychties_outgroups.prot XP_038672178.1_Scyliorhinus-canicula-tmt3 >> tmt3.prot
samtools faidx ../Chondrychties_outgroups.prot XP_051874554.1_Pristis-pectinat-tmt3 >> tmt3.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align tmt3.prot -output tmt3.aln
xargs samtools faidx tmt3.aln < tmt3 > tmt3.fish.aln
trimal -in tmt3.fish.aln -backtrans tmt3.fa -automated1 -out tmt3.fish.cds.trimmed.aln
trimal -in tmt3.aln -automated1 -out tmt3.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh tmt3.prot.trimmed.aln
echo "XP_038672178.1_Scyliorhinus-canicula-tmt3" > outgroup.id
echo "XP_051874554.1_Pristis-pectinat-tmt3" >> outgroup.id
xargs samtools faidx tmt3.prot.trimmed.aln < tmt3 > tmt3.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R tmt3.prot.trimmed.aln.treefile outgroup.id tmt3.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R tmt3.prot.trimmed.aln.nooutgroup.treefile tmt3.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R tmt3.prot.trimmed.aln.nooutgroup.treefile tmt3.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## Va ######## 


mkdir va ; cd va
cp ../Clades/va ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < va > va.prot
xargs samtools faidx ../Complete_opsins.fa < va > va.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007894489.1_Callorhinchus-milii-va >> va.prot
samtools faidx ../Chondrychties_outgroups.prot XP_041065233.1_Carcharodon-carcharias-va >> va.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align va.prot -output va.aln
xargs samtools faidx va.aln < va > va.fish.aln
trimal -in va.fish.aln -backtrans va.fa -automated1 -out va.fish.cds.trimmed.aln
trimal -in va.aln -automated1 -out va.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh va.prot.trimmed.aln
echo "XP_007894489.1_Callorhinchus-milii-va" > outgroup.id
echo "XP_041065233.1_Carcharodon-carcharias-va" >> outgroup.id
xargs samtools faidx va.prot.trimmed.aln < va > va.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R va.prot.trimmed.aln.treefile outgroup.id va.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R va.prot.trimmed.aln.nooutgroup.treefile va.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R va.prot.trimmed.aln.nooutgroup.treefile va.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## Opn9 ########


mkdir opn9 ; cd opn9
cp ../Clades/opn9 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn9 > opn9.prot
xargs samtools faidx ../Complete_opsins.fa < opn9 > opn9.fa
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn9.prot -output opn9.aln

samtools faidx ../Chondrychties_outgroups.prot XP_007901513.1_Callorhinchus_milii-opn5 > shark_opn9.prot

mafft --add shark_opn9.prot --keeplength opn9.aln > temp.aln ; mv temp.aln opn9.aln

xargs samtools faidx opn9.aln < opn9 > opn9.fish.aln
trimal -in opn9.fish.aln -backtrans opn9.fa -automated1 -out opn9.fish.cds.trimmed.aln
trimal -in opn9.aln -automated1 -out opn9.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn9.prot.trimmed.aln
echo "XP_007901513.1_Callorhinchus_milii-opn5" > outgroup.id
xargs samtools faidx opn9.prot.trimmed.aln < opn9 > opn9.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn9.prot.trimmed.aln.treefile outgroup.id opn9.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn9.prot.trimmed.aln.nooutgroup.treefile opn9.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn9.prot.trimmed.aln.nooutgroup.treefile opn9.prot.trimmed.aln.nooutgroup.collapsed95.treefile

########## Opn5 ########


mkdir opn5 ; cd opn5
cp ../Clades/opn5 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn5 > opn5.prot
xargs samtools faidx ../Complete_opsins.fa < opn5 > opn5.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007901513.1_Callorhinchus_milii-opn5 >> opn5.prot
samtools faidx ../Chondrychties_outgroups.prot XP_048453349.1_Rhincodon_typus-opn5 >> opn5.prot


/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn5.prot -output opn5.aln
xargs samtools faidx opn5.aln < opn5 > opn5.fish.aln
trimal -in opn5.fish.aln -backtrans opn5.fa -automated1 -out opn5.fish.cds.trimmed.aln
trimal -in opn5.aln -automated1 -out opn5.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn5.prot.trimmed.aln
echo "XP_007901513.1_Callorhinchus_milii-opn5" > outgroup.id
echo "XP_048453349.1_Rhincodon_typus-opn5" >> outgroup.id
xargs samtools faidx opn5.prot.trimmed.aln < opn5 > opn5.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn5.prot.trimmed.aln.treefile outgroup.id opn5.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn5.prot.trimmed.aln.nooutgroup.treefile opn5.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn5.prot.trimmed.aln.nooutgroup.treefile opn5.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## Opn6 ######## 


mkdir opn6 ; cd opn6
cp ../Clades/opn6 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn6 > opn6.prot
xargs samtools faidx ../Complete_opsins.fa < opn6 > opn6.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007887390.1_Callorhinchus_milii-opn6 >> opn6.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn6.prot -output opn6.aln
xargs samtools faidx opn6.aln < opn6 > opn6.fish.aln
trimal -in opn6.fish.aln -backtrans opn6.fa -automated1 -out opn6.fish.cds.trimmed.aln
trimal -in opn6.aln -automated1 -out opn6.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn6.prot.trimmed.aln
echo "XP_007887390.1_Callorhinchus_milii-opn6" > outgroup.id
xargs samtools faidx opn6.prot.trimmed.aln < opn6 > opn6.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn6.prot.trimmed.aln.treefile outgroup.id opn6.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn6.prot.trimmed.aln.nooutgroup.treefile opn6.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn6.prot.trimmed.aln.nooutgroup.treefile opn6.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## Opn7a ######## 


mkdir opn7a ; cd opn7a
cp ../Clades/opn7a ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn7a > opn7a.prot
xargs samtools faidx ../Complete_opsins.fa < opn7a > opn7a.fa


shuf -n10 ../../Phylogenies_Opsins/ALL_opsins_tree/OLD_CLADES/opn7a > random_id
sort random_id > temp ; mv temp random_id
grep ">" opn7a.fa | sed 's/>//g' | sort > all.id
comm -23 all.id random_id > rest_of_sequences
xargs samtools faidx opn7a.prot < random_id > random.prot
xargs samtools faidx opn7a.prot < rest_of_sequences > rest_of_sequences.prot
samtools faidx ../Chondrychties_outgroups.prot XP_007883754.2_Callorhinchus_milii-opn7a >> rest_of_sequences.prot


/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align random.prot -output random.aln

mafft --add rest_of_sequences.prot --keeplength random.aln > opn7a.aln

xargs samtools faidx opn7a.aln < opn7a > opn7a.fish.aln
trimal -in opn7a.fish.aln -backtrans opn7a.fa -automated1 -out opn7a.fish.cds.trimmed.aln
trimal -in opn7a.aln -automated1 -out opn7a.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn7a.prot.trimmed.aln
echo "XP_007883754.2_Callorhinchus_milii-opn7a" > outgroup.id
xargs samtools faidx opn7a.prot.trimmed.aln < opn7a > opn7a.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn7a.prot.trimmed.aln.treefile outgroup.id opn7a.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn7a.prot.trimmed.aln.nooutgroup.treefile opn7a.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn7a.prot.trimmed.aln.nooutgroup.treefile opn7a.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## Opn7b ######## 


mkdir opn7b ; cd opn7b
cp ../Clades/opn7b ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn7b > opn7b.prot
xargs samtools faidx ../Complete_opsins.fa < opn7b > opn7b.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007883754.2_Callorhinchus_milii-opn7a >> opn7b.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn7b.prot -output opn7b.aln
xargs samtools faidx opn7b.aln < opn7b > opn7b.fish.aln
trimal -in opn7b.fish.aln -backtrans opn7b.fa -automated1 -out opn7b.fish.cds.trimmed.aln
trimal -in opn7b.aln -automated1 -out opn7b.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn7b.prot.trimmed.aln
echo "XP_007883754.2_Callorhinchus_milii-opn7a" > outgroup.id
xargs samtools faidx opn7b.prot.trimmed.aln < opn7b > opn7b.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn7b.prot.trimmed.aln.treefile outgroup.id opn7b.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn7b.prot.trimmed.aln.nooutgroup.treefile opn7b.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn7b.prot.trimmed.aln.nooutgroup.treefile opn7b.prot.trimmed.aln.nooutgroup.collapsed95.treefile


########## opn4m1_3 ######## 


mkdir opn4m1_3 ; cd opn4m1_3
cp ../Clades/opn4m1_3 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn4m1_3 > opn4m1_3.prot
xargs samtools faidx ../Complete_opsins.fa < opn4m1_3 > opn4m1_3.fa
samtools faidx ../Chondrychties_outgroups.prot NP_001279357.1_Callorhinchus-milii-opn4m1_3 >> opn4m1_3.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn4m1_3.prot -output opn4m1_3.aln
xargs samtools faidx opn4m1_3.aln < opn4m1_3 > opn4m1_3.fish.aln
trimal -in opn4m1_3.fish.aln -backtrans opn4m1_3.fa -automated1 -out opn4m1_3.fish.cds.trimmed.aln
trimal -in opn4m1_3.aln -automated1 -out opn4m1_3.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn4m1_3.prot.trimmed.aln
echo "NP_001279357.1_Callorhinchus-milii-opn4m1_3" > outgroup.id
xargs samtools faidx opn4m1_3.prot.trimmed.aln < opn4m1_3 > opn4m1_3.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn4m1_3.prot.trimmed.aln.treefile outgroup.id opn4m1_3.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn4m1_3.prot.trimmed.aln.nooutgroup.treefile opn4m1_3.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn4m1_3.prot.trimmed.aln.nooutgroup.treefile opn4m1_3.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## lws ######## 


mkdir lws ; cd lws
cp ../Clades/lws ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < lws > lws.prot
xargs samtools faidx ../Complete_opsins.fa < lws > lws.fa


grep "Polypterus\|calabaricus\|Lepiso\|calva\|Danio_rerio\|Astyanax\|Megalops\|Oryzias" ../../Phylogenies_Opsins/ALL_opsins_tree/OLD_CLADES/lws > random_id
sort random_id > temp ; mv temp random_id
grep ">" lws.fa | sed 's/>//g' | sort > all.id
comm -23 all.id random_id > rest_of_sequences
xargs samtools faidx lws.prot < random_id > random.prot
xargs samtools faidx lws.prot < rest_of_sequences > rest_of_sequences.prot
cat ../outgroup_lws.prot >> random.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align random.prot -output random.aln

mafft --add rest_of_sequences.prot --keeplength random.aln > lws.aln

xargs samtools faidx lws.aln < lws > lws.fish.aln
trimal -in lws.fish.aln -backtrans lws.fa -automated1 -out lws.fish.cds.trimmed.aln

cp lws.aln lws.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh lws.prot.trimmed.aln
grep ">"  ../outgroup_lws.prot | sed 's/>//g' > outgroup.id
xargs samtools faidx lws.prot.trimmed.aln < lws > lws.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R lws.prot.trimmed.aln.treefile outgroup.id lws.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R lws.prot.trimmed.aln.nooutgroup.treefile lws.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R lws.prot.trimmed.aln.nooutgroup.treefile lws.prot.trimmed.aln.nooutgroup.collapsed95.treefile


grep -v "Lethenteron_camtschaticum" outgroup.id > toremove.id
Rscript Remove_outgroup.R lws.prot.trimmed.aln.treefile toremove.id lws.prot.trimmed.aln.treefile.BETTER


cp lws.prot.trimmed.aln.treefile tonotuse_lws.prot.trimmed.aln.treefile
mv lws.prot.trimmed.aln.treefile.BETTER lws.prot.trimmed.aln.treefile



########## rgr ######## 


mkdir rgr ; cd rgr
cp ../Clades/rgr ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < rgr > rgr.prot
xargs samtools faidx ../Complete_opsins.fa < rgr > rgr.fa
samtools faidx ../Chondrychties_outgroups.prot XM_007898380.1_Callorhinchus-milii-rgr >> rgr.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align rgr.prot -output rgr.aln
xargs samtools faidx rgr.aln < rgr > rgr.fish.aln
trimal -in rgr.fish.aln -backtrans rgr.fa -automated1 -out rgr.fish.cds.trimmed.aln
trimal -in rgr.aln -automated1 -out rgr.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh rgr.prot.trimmed.aln
echo "XM_007898380.1_Callorhinchus-milii-rgr" > outgroup.id
xargs samtools faidx rgr.prot.trimmed.aln < rgr > rgr.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R rgr.prot.trimmed.aln.treefile outgroup.id rgr.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R rgr.prot.trimmed.aln.nooutgroup.treefile rgr.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R rgr.prot.trimmed.aln.nooutgroup.treefile rgr.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## rh1 ######## 


mkdir rh1 ; cd rh1
cp ../Clades/rh1 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < rh1 > rh1.prot
xargs samtools faidx ../Complete_opsins.fa < rh1 > rh1.fa
samtools faidx ../Chondrychties_outgroups.prot A0A7T8R2L6_Rhincodon-typus-rh1 >> rh1.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align rh1.prot -output rh1.aln
xargs samtools faidx rh1.aln < rh1 > rh1.fish.aln
trimal -in rh1.fish.aln -backtrans rh1.fa -automated1 -out rh1.fish.cds.trimmed.aln
trimal -in rh1.aln -automated1 -out rh1.prot.trimmed.aln
sbatch --qos=6hours -c 20 --mem=50G iqtree_phylo.sh rh1.prot.trimmed.aln
echo "A0A7T8R2L6_Rhincodon-typus-rh1" > outgroup.id
xargs samtools faidx rh1.prot.trimmed.aln < rh1 > rh1.fish.prot.trimmed.aln

cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R rh1.prot.trimmed.aln.treefile outgroup.id rh1.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R rh1.prot.trimmed.aln.nooutgroup.treefile rh1.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R rh1.prot.trimmed.aln.nooutgroup.treefile rh1.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## rh2 ######## 
 

mkdir rh2 ; cd rh2
cp ../Clades/rh2 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < rh2 > rh2.prot
xargs samtools faidx ../Complete_opsins.fa < rh2 > rh2.fa
samtools faidx ../Chondrychties_outgroups.prot XP_051890468.1_Pristis-pectinata--rh2 >> rh2.prot
samtools faidx ../Chondrychties_outgroups.prot XP_059806431.1_Hypanus-sabinus-rh2 >> rh2.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align rh2.prot -output rh2.aln
xargs samtools faidx rh2.aln < rh2 > rh2.fish.aln
trimal -in rh2.fish.aln -backtrans rh2.fa -automated1 -out rh2.fish.cds.trimmed.aln
trimal -in rh2.aln -automated1 -out rh2.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh rh2.prot.trimmed.aln
echo "XP_051890468.1_Pristis-pectinata--rh2" > outgroup.id
echo "XP_059806431.1_Hypanus-sabinus-rh2" >> outgroup.id
xargs samtools faidx rh2.prot.trimmed.aln < rh2 > rh2.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R rh2.prot.trimmed.aln.treefile outgroup.id rh2.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R rh2.prot.trimmed.aln.nooutgroup.treefile rh2.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R rh2.prot.trimmed.aln.nooutgroup.treefile rh2.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## rrh ######## 


mkdir rrh ; cd rrh
cp ../Clades/rrh ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < rrh > rrh.prot
xargs samtools faidx ../Complete_opsins.fa < rrh > rrh.fa
samtools faidx ../Chondrychties_outgroups.prot XP_051866389.1_Pristis-pectinata--rrh >> rrh.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align rrh.prot -output rrh.aln
xargs samtools faidx rrh.aln < rrh > rrh.fish.aln
trimal -in rrh.fish.aln -backtrans rrh.fa -automated1 -out rrh.fish.cds.trimmed.aln
trimal -in rrh.aln -automated1 -out rrh.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh rrh.prot.trimmed.aln
echo "XP_051866389.1_Pristis-pectinata--rrh" > outgroup.id
xargs samtools faidx rrh.prot.trimmed.aln < rrh > rrh.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R rrh.prot.trimmed.aln.treefile outgroup.id rrh.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R rrh.prot.trimmed.aln.nooutgroup.treefile rrh.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R rrh.prot.trimmed.aln.nooutgroup.treefile rrh.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## opn3 ######## 


mkdir opn3 ; cd opn3
cp ../Clades/opn3 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn3 > opn3.prot
xargs samtools faidx ../Complete_opsins.fa < opn3 > opn3.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007892106.1_Callorhinchus-milii-opn3 >> opn3.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn3.prot -output opn3.aln
xargs samtools faidx opn3.aln < opn3 > opn3.fish.aln
trimal -in opn3.fish.aln -backtrans opn3.fa -automated1 -out opn3.fish.cds.trimmed.aln
trimal -in opn3.aln -automated1 -out opn3.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn3.prot.trimmed.aln
echo "XP_007892106.1_Callorhinchus-milii-opn3" > outgroup.id
xargs samtools faidx opn3.prot.trimmed.aln < opn3 > opn3.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn3.prot.trimmed.aln.treefile outgroup.id opn3.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn3.prot.trimmed.aln.nooutgroup.treefile opn3.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn3.prot.trimmed.aln.nooutgroup.treefile opn3.prot.trimmed.aln.nooutgroup.collapsed95.treefile



########## sws2 ######## ==> (mtInv+F+R8) 

cp ../Clades/sws2 ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < sws2 > sws2.prot
xargs samtools faidx ../Complete_opsins.fa < sws2 > sws2.fa


grep "Polypterus\|calabaricus\|Lepiso\|calva\|Danio_rerio\|Astyanax\|Megalops\|Oryzias" ../../Phylogenies_Opsins/ALL_opsins_tree/OLD_CLADES/sws2 > random_id
sort random_id > temp ; mv temp random_id
grep ">" sws2.fa | sed 's/>//g' | sort > all.id
comm -23 all.id random_id > rest_of_sequences
xargs samtools faidx sws2.prot < random_id > random.prot
xargs samtools faidx sws2.prot < rest_of_sequences > rest_of_sequences.prot
cat ../sws1_lamprey.prot >> rest_of_sequences.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align random.prot -output random.aln

mafft --add rest_of_sequences.prot --keeplength random.aln > sws2.aln

xargs samtools faidx sws2.aln < sws2 > sws2.fish.aln
trimal -in sws2.fish.aln -backtrans sws2.fa -automated1 -out sws2.fish.cds.trimmed.aln
trimal -in sws2.aln -automated1 -out sws2.prot.trimmed.aln


#sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh sws2.aln
#sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh sws2.prot.trimmed.aln

echo "Geotria_australis_AY366495_sws1" > outgroup.id
xargs samtools faidx sws2.prot.trimmed.aln < sws2 > sws2.fish.prot.trimmed.aln

sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh sws2.fish.prot.trimmed.aln

cp sws2.fish.prot.trimmed.aln.treefile sws2.prot.trimmed.aln.treefile

#cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
#Rscript Remove_outgroup.R sws2.prot.trimmed.aln.treefile outgroup.id sws2.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R sws2.prot.trimmed.aln.treefile sws2.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R sws2.prot.trimmed.aln.treefile sws2.prot.trimmed.aln.nooutgroup.collapsed95.treefile
cp sws2.prot.trimmed.aln.treefile sws2.prot.trimmed.aln.nooutgroup.treefile




########## opn8a ######## =


mkdir opn8a ; cd opn8a
cp ../Clades/opn8a ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn8a > opn8a.prot
xargs samtools faidx ../Complete_opsins.fa < opn8a > opn8a.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007901512.2_Callorhinchus_milii-opn8b >> opn8a.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn8a.prot -output opn8a.aln
xargs samtools faidx opn8a.aln < opn8a > opn8a.fish.aln
trimal -in opn8a.fish.aln -backtrans opn8a.fa -automated1 -out opn8a.fish.cds.trimmed.aln
trimal -in opn8a.aln -automated1 -out opn8a.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn8a.prot.trimmed.aln
echo "XP_007901512.2_Callorhinchus_milii-opn8b" > outgroup.id
xargs samtools faidx opn8a.prot.trimmed.aln < opn8a > opn8a.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn8a.prot.trimmed.aln.treefile outgroup.id opn8a.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn8a.prot.trimmed.aln.nooutgroup.treefile opn8a.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn8a.prot.trimmed.aln.nooutgroup.treefile opn8a.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## opn8b ######## =


mkdir opn8b ; cd opn8b
cp ../Clades/opn8b ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn8b > opn8b.prot
xargs samtools faidx ../Complete_opsins.fa < opn8b > opn8b.fa
samtools faidx ../Chondrychties_outgroups.prot XP_007901512.2_Callorhinchus_milii-opn8b >> opn8b.prot
/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align opn8b.prot -output opn8b.aln
xargs samtools faidx opn8b.aln < opn8b > opn8b.fish.aln
trimal -in opn8b.fish.aln -backtrans opn8b.fa -automated1 -out opn8b.fish.cds.trimmed.aln
trimal -in opn8b.aln -automated1 -out opn8b.prot.trimmed.aln
sbatch --qos=1day -c 20 --mem=50G iqtree_phylo.sh opn8b.prot.trimmed.aln
echo "XP_007901512.2_Callorhinchus_milii-opn8b" > outgroup.id
xargs samtools faidx opn8b.prot.trimmed.aln < opn8b > opn8b.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
Rscript Remove_outgroup.R opn8b.prot.trimmed.aln.treefile outgroup.id opn8b.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn8b.prot.trimmed.aln.nooutgroup.treefile opn8b.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn8b.prot.trimmed.aln.nooutgroup.treefile opn8b.prot.trimmed.aln.nooutgroup.collapsed95.treefile




########## opn8c ######## 


mkdir opn8c ; cd opn8c
cp ../Clades/opn8c ./
cp ../iqtree_phylo.sh ./
xargs samtools faidx ../Complete_opsins.prot < opn8c > opn8c.prot
xargs samtools faidx ../Complete_opsins.fa < opn8c > opn8c.fa


grep "Polypterus\|calabaricus\|Lepiso\|calva\|Danio_rerio" ../../Phylogenies_Opsins/ALL_opsins_tree/OLD_CLADES/opn8c > random_id
sort random_id > temp ; mv temp random_id
grep ">" opn8c.fa | sed 's/>//g' | sort > all.id
comm -23 all.id random_id > rest_of_sequences
xargs samtools faidx opn8c.prot < random_id > random.prot
xargs samtools faidx opn8c.prot < rest_of_sequences > rest_of_sequences.prot


/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align random.prot -output random.aln

mafft --add rest_of_sequences.prot --keeplength random.aln > opn8c.aln

xargs samtools faidx opn8c.aln < opn8c > opn8c.fish.aln
trimal -in opn8c.fish.aln -backtrans opn8c.fa -automated1 -out opn8c.fish.cds.trimmed.aln
cp  opn8c.aln opn8c.prot.trimmed.aln
sbatch --qos=6hours -c 20 --mem=50G iqtree_phylo.sh opn8c.prot.trimmed.aln
xargs samtools faidx opn8c.prot.trimmed.aln < opn8c > opn8c.fish.prot.trimmed.aln


cp ../Remove_outgroup.R ./ 
cp ../Collapse_tree_90.R ./ 
cp ../Collapse_tree_95.R ./ 
cp opn8c.prot.trimmed.aln.treefile opn8c.prot.trimmed.aln.nooutgroup.treefile
Rscript Collapse_tree_90.R opn8c.prot.trimmed.aln.treefile opn8c.prot.trimmed.aln.nooutgroup.collapsed90.treefile
Rscript Collapse_tree_95.R opn8c.prot.trimmed.aln.treefile opn8c.prot.trimmed.aln.nooutgroup.collapsed95.treefile







##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Make a nice table of the nb of gene per species, as well as the nb of pseudo, incomplete ... ######################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#BUSCO90teleost_BUSCO80nonteleost_species.txt => file with the lest of species in the dataset ...

cd /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/
cp ../BUSCO_results/BUSCO90teleost_BUSCO80nonteleost_species.txt ./
sed -i 's/Taenioides_sp._WSHXM2023/Taenioides_sp_WSHXM2023/g' BUSCO90teleost_BUSCO80nonteleost_species.txt

rm ALL_opsins_classifications.csv

for species in `cat BUSCO90teleost_BUSCO80nonteleost_species.txt` ; do 

	grep ">$species---" Complete_opsins.fa  | sed 's/>//g' > current_sp_complete.id

	for gene in `cat current_sp_complete.id` ; do 

		clade=`grep "$gene" Clades/* | sed 's/:.*//g' | sed 's/.*\///g'`


		echo "Complete,$gene,$clade,$species" >> ALL_opsins_classifications.csv


	done 

done


for species in `cat BUSCO90teleost_BUSCO80nonteleost_species.txt` ; do 

	grep ">$species---" Pseudogene_opsins.fa  | sed 's/>//g' > current_sp_pseudo.id
	
	for gene in `cat current_sp_pseudo.id` ; do 

		samtools faidx Pseudogene_opsins.fa $gene > current_pseudo.fa 

		blastx -query current_pseudo.fa -db Complete_opsins.prot -outfmt "6 sseqid" -num_threads 8 -max_target_seqs 1 -out blastx_result
		best_blastx=`head -1 blastx_result`

		clade=`grep "$best_blastx" ALL_opsins_classifications.csv | cut -f3 -d ","`
		
		echo "Pseudogene,$gene,$clade,$species" >> ALL_opsins_classifications.csv

	done 


done


for species in `cat BUSCO90teleost_BUSCO80nonteleost_species.txt` ; do 

	grep ">$species---" Incomplete_opsins.fa  | sed 's/>//g' > current_sp_pseudo.id
	
	for gene in `cat current_sp_pseudo.id` ; do 

		samtools faidx Incomplete_opsins.fa $gene > current_pseudo.fa 

		blastx -query current_pseudo.fa -db Complete_opsins.prot -outfmt "6 sseqid" -num_threads 8 -max_target_seqs 1 -out blastx_result
		best_blastx=`head -1 blastx_result`

		clade=`grep "$best_blastx" ALL_opsins_classifications.csv | cut -f3 -d ","`
		
		echo "Incomplete,$gene,$clade,$species" >> ALL_opsins_classifications.csv

	done 


done


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract some infos ################################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


rm Subfamilies_stats.csv 
for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`


	echo "$clade_name" 


	xargs samtools faidx Complete_opsins.prot < Clades/$clade_name > current.prot
	rm current.prot.fai ; samtools faidx current.prot

	mean_length=`awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' current.prot.fai`
	mean_exon_nb=`sed 's/.*---//g' Clades/$clade_name | sed 's/_exons//g' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`

	echo "$clade_name,$mean_length,$mean_exon_nb" >> Subfamilies_stats.csv 

done


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Species tree / Gene tree reconciliation ###########################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#I will use two tree : one made with IQ-tree and one made with ASTRAL (See script MakeTree_Actino..). 

cd /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/

conda activate treerecs_env


rm -r ASTRAL_tree_reconcilation_95 ; rm -r IQTREE_tree_reconcilation_95
rm -r ASTRAL_tree_reconcilation_90 ; rm -r IQTREE_tree_reconcilation_90
rm -r ASTRAL_tree_reconcilation_nocollapse ; rm -r IQTREE_tree_reconcilation_nocollapse



for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`

	echo "$clade_name" 

	
	treerecs -s /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/ASTRAL/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel -g $clade_name/$clade_name.prot.trimmed.aln.nooutgroup.treefile ‐‐reroot -r -v -f ‐‐fevent --outdir ASTRAL_tree_reconcilation_nocollapse
	treerecs -s /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/Concatenated_Alignments/LGmodel/Dating_min_01/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel -g $clade_name/$clade_name.prot.trimmed.aln.nooutgroup.treefile ‐‐reroot -r -v -f ‐‐fevent --outdir IQTREE_tree_reconcilation_nocollapse


done > Treerecs_error_msg_nocollapse



for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`

	echo "$clade_name" 

	
	treerecs -s /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/ASTRAL/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel -g $clade_name/$clade_name.prot.trimmed.aln.nooutgroup.collapsed90.treefile ‐‐reroot -r -v -f ‐‐fevent --outdir ASTRAL_tree_reconcilation_90
	treerecs -s /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/Concatenated_Alignments/LGmodel/Dating_min_01/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel -g $clade_name/$clade_name.prot.trimmed.aln.nooutgroup.collapsed90.treefile ‐‐reroot -r -v -f ‐‐fevent --outdir IQTREE_tree_reconcilation_90


done > Treerecs_error_msg_90



for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`

	echo "$clade_name" 

	
	treerecs -s /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/ASTRAL/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel -g $clade_name/$clade_name.prot.trimmed.aln.nooutgroup.collapsed95.treefile ‐‐reroot -r -v -f ‐‐fevent --outdir ASTRAL_tree_reconcilation_95
	treerecs -s /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/Concatenated_Alignments/LGmodel/Dating_min_01/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel -g $clade_name/$clade_name.prot.trimmed.aln.nooutgroup.collapsed95.treefile ‐‐reroot -r -v -f ‐‐fevent --outdir IQTREE_tree_reconcilation_95


done > Treerecs_error_msg_95





#Reformat treerecs output and compute birth and death events with NOTUNG --- ASTRAL TREE


for folder in ASTRAL_tree_reconcilation* ; do 


	folder_name=`echo "$folder" | sed 's/\///g'`
	cd $folder 

	for i in *.nwk ; do grep -v ">" $i > $i.TREE  ; done

	rm -r TreeRecs_Summary ; mkdir TreeRecs_Summary
	cd TreeRecs_Summary
	cp ../*.TREE ./
	cp /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/ASTRAL/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel ./

	sed -i 's/Carassius_auratus_red_var/Carassius_red_var/g' *.TREE
	sed -i 's/Carassius_auratus_red_var/Carassius_red_var/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel
	sed -i 's/Cyprinus_carpio_carpio/Cyprinus_second/g' *.TREE
	sed -i 's/Cyprinus_carpio_carpio/Cyprinus_second/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel


	echo "/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/$folder_name/TreeRecs_Summary/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel"  > batch.txt
	ls -l | grep ".TREE$" | sed 's/.* //g' >> batch.txt

	java -jar /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/GeneTree_SpeciesTree_Reconcilitation/NOTUNG/Notung-2.9.1.5.jar -b batch.txt --reconcile --speciestag prefix --phylogenomics 

	sed -i 's/.prot.aln.trimal.nooutgroup.collapsed.treefile_recs.nwk.TREE//g' *.txt

	sed -i 's/Carassius_red_var/Carassius_auratus_red_var/g' *.txt
	sed -i 's/Cyprinus_second/Cyprinus_carpio_carpio/g' *.txt

	sed -i 's/Carassius_red_var/Carassius_auratus_red_var/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel
	sed -i 's/Cyprinus_second/Cyprinus_carpio_carpio/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel


	cat /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/$folder_name/TreeRecs_Summary/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel | tr , '\n' | tr ')' '\n' | sed 's/(//g' | sed 's/:/ /g' | grep -v 'Node1;' > Species_branchlength.tsv

	cd ../../

done




for folder in IQTREE_tree_reconcilation* ; do 


	folder_name=`echo "$folder" | sed 's/\///g'`
	cd $folder 

	for i in *.nwk ; do grep -v ">" $i > $i.TREE  ; done

	rm -r TreeRecs_Summary ; mkdir TreeRecs_Summary
	cd TreeRecs_Summary
	cp ../*.TREE ./
	cp /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Species_Phylogeny/Concatenated_Alignments/LGmodel/Dating_min_01/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel ./

	sed -i 's/Carassius_auratus_red_var/Carassius_red_var/g' *.TREE
	sed -i 's/Carassius_auratus_red_var/Carassius_red_var/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel
	sed -i 's/Cyprinus_carpio_carpio/Cyprinus_second/g' *.TREE
	sed -i 's/Cyprinus_carpio_carpio/Cyprinus_second/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel


	echo "/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/$folder_name/TreeRecs_Summary/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel"  > batch.txt
	ls -l | grep ".TREE$" | sed 's/.* //g' >> batch.txt

	java -jar /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/GeneTree_SpeciesTree_Reconcilitation/NOTUNG/Notung-2.9.1.5.jar -b batch.txt --reconcile --speciestag prefix --phylogenomics 

	sed -i 's/.prot.aln.trimal.nooutgroup.collapsed.treefile_recs.nwk.TREE//g' *.txt

	sed -i 's/Carassius_red_var/Carassius_auratus_red_var/g' *.txt
	sed -i 's/Cyprinus_second/Cyprinus_carpio_carpio/g' *.txt

	sed -i 's/Carassius_red_var/Carassius_auratus_red_var/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel
	sed -i 's/Cyprinus_second/Cyprinus_carpio_carpio/g' AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel


	cat /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/$folder_name/TreeRecs_Summary/AMAS_concatenated_alignment_2000BUSCO.fa.timetree.nwk.nodelabel | tr , '\n' | tr ')' '\n' | sed 's/(//g' | sed 's/:/ /g' | grep -v 'Node1;' > Species_branchlength.tsv

	cd ../../

done



##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Codeml -- Dup vs Spec ##############################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================




cd Phylogenies_Opsins_v2/
mkdir Reconciled_trees_95 ; cd Reconciled_trees_95


cp /scicore/home/salzburg/polica0000/Vertebrates_Taste_Receptors/Birds_dNdS_V1R/OneRatio.ctl ./
cp /scicore/home/salzburg/polica0000/Vertebrates_Taste_Receptors/Birds_dNdS_V1R/TwoRatio.ctl ./
cp ../prepare_codeml_dup.R ./

for clade in ../Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/..\/Clades\///g'` 


	cp ../ASTRAL_tree_reconcilation_95/TreeRecs_Summary/$clade_name.prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE.reconciled ./

	head -1 $clade_name.prot.trimmed.aln.nooutgroup.collapsed95.treefile_recs.nwk.TREE.reconciled > $clade_name.95.reconciled

	Rscript prepare_codeml_dup.R $clade_name.95.reconciled $clade_name.dup_labels.txt $clade_name.all_tips_labels $clade_name.all_node_labels $clade_name.myphylo_noloss.nwk

	nb_seq=`wc -l < $clade_name.all_tips_labels`
	
	for i in `seq 1 $nb_seq` ; do
		echo "sequence_$i" >> $clade_name.new_labels
	done
	
	paste -d "," $clade_name.all_tips_labels $clade_name.new_labels > $clade_name.renaming_file.csv
	
	
	cp ../$clade_name/$clade_name.fish.cds.trimmed.aln ./
	sed -i 's/Carassius_auratus_red_var/Carassius_red_var/g' $clade_name.fish.cds.trimmed.aln
	sed -i 's/Cyprinus_carpio_carpio/Cyprinus_second/g' $clade_name.fish.cds.trimmed.aln
	
	for line in `cat $clade_name.renaming_file.csv` ; do 
		old_name=`echo "$line" | cut -f1 -d ","`
		new_name=`echo "$line" | cut -f2 -d ","`
	
	
		sed -i "s/$old_name/$new_name/g" $clade_name.myphylo_noloss.nwk
		sed -i "s/$old_name/$new_name/g" $clade_name.fish.cds.trimmed.aln
	
	done

	tac $clade_name.all_node_labels > temp.txt
	mv temp.txt $clade_name.all_node_labels


	for nodes in `cat $clade_name.all_node_labels` ; do sed -i "s/$nodes//g" $clade_name.myphylo_noloss.nwk ; done

	sed -i 's/_#1/ #1/g' $clade_name.myphylo_noloss.nwk


	mkdir $clade_name ; mkdir $clade_name/TwoRatio ; mkdir $clade_name/OneRatio
	cp $clade_name.fish.cds.trimmed.aln $clade_name/OneRatio/ ; cp $clade_name.myphylo_noloss.nwk $clade_name/OneRatio/
	cp $clade_name.fish.cds.trimmed.aln $clade_name/TwoRatio/ ; cp $clade_name.myphylo_noloss.nwk $clade_name/TwoRatio/

	sed "s/mysequences/$clade_name.fish.cds.trimmed.aln/g" TwoRatio.ctl | sed "s/mytreefile/$clade_name.myphylo_noloss.nwk/g" | sed 's/myoutput.mcl/TwoRatio.mcl/g' > $clade_name/TwoRatio/TwoRatio.ctl
	sed "s/mysequences/$clade_name.fish.cds.trimmed.aln/g" OneRatio.ctl | sed "s/mytreefile/$clade_name.myphylo_noloss.nwk/g" | sed 's/myoutput.mcl/OneRatio.mcl/g' > $clade_name/OneRatio/OneRatio.ctl


	cd $clade_name/OneRatio 
	sbatch --job-name=PAML_ORM --qos=2weeks -c 3 --mem=5G --wrap='eval "$(conda shell.bash hook)" ; conda activate PAML ; codeml OneRatio.ctl'

	cd ../TwoRatio
	sbatch --job-name=PAML_TRM --qos=2weeks -c 3 --mem=5G --wrap='eval "$(conda shell.bash hook)" ; conda activate PAML ; codeml TwoRatio.ctl'


	cd ../../

	cp $clade_name* $clade_name/

	rm $clade_name*

done 





### Now do the same but with BUSTED


for clade in ../Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/..\/Clades\///g'` 


	cd $clade_name 
	rm -r BUSTED ; mkdir BUSTED
	cd BUSTED

	cp ../$clade_name.fish.cds.trimmed.aln ./
	cp ../$clade_name.myphylo_noloss.nwk ./
	cp ../../run_busted.sh ./

	sed -i 's/ #1/{TestBranch}/g' $clade_name.myphylo_noloss.nwk


	sbatch -c 20 --qos=6hours --mem=30G --job-name=BUSTED.$clade_name -e error.busted.$clade_name.out -o slurm.busted.$clade_name.out run_busted.sh $clade_name

	cd ../../

done

mkdir BUSTED 

for clade in ../Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/..\/Clades\///g'` 

	cp $clade_name/BUSTED/*.json BUSTED/

done






### Codeml takes way too much time, so one might prefer to just run BUSTED


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Compute the omega ratio per branch for each gene tree #############################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


conda activate HyPhy_Env





#### Using the Gene trees
for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`

	echo "$clade_name" 

	cd $clade_name

	cp ../run_fitmg94.sh ./ 
	cp ../FitMG94.bf ./


	sbatch -c 10 --qos=1week --mem=20G --job-name=fitMG -e error.fitMG.$clade_name.out -o slurm.gitMG.$clade_name.out run_fitmg94.sh $clade_name

	cd ../


done



#### Using the reconciled gene trees

cd Phylogenies_Opsins_v2/Reconciled_trees_95/

for clade in ../Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/..\/Clades\///g'` 


	cd $clade_name 
	rm -r Omega_per_branch ; mkdir Omega_per_branch
	cd Omega_per_branch

	cp ../$clade_name.fish.cds.trimmed.aln ./
	cp ../$clade_name.myphylo_noloss.nwk ./
	cp ../../../run_fitmg94.sh ./
	cp ../../../FitMG94.bf ./

	sed 's/ #1//g' $clade_name.myphylo_noloss.nwk > $clade_name.prot.trimmed.aln.nooutgroup.treefile

	sbatch -c 10 --qos=1week --mem=15G --job-name=fitMG94.$clade_name -e error.fitmg94.$clade_name.out -o slurm.fitmg94.$clade_name.out run_fitmg94.sh $clade_name

	cd ../../

done




for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`
	echo "$clade_name"
	wc -l $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json

done



# Extract Results to tables ! 


cd /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2

rm -r OMEGA_per_branch_genetree ; mkdir OMEGA_per_branch_genetree
rm -r OMEGA_per_branch_speciestree ; mkdir OMEGA_per_branch_speciestree


for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`

	grep "LB\":" $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
	grep "MLE\":" $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
	grep "UB\":" $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
	grep "\"dN\"" $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
	grep "\"dS\"" $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
	grep -B2 "LB\":" $clade_name/$clade_name.fish.cds.trimmed.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels
	sed -i "s/^/$clade_name,/g" curr_labels


	paste -d "," curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_labels curr_dN_values.txt curr_dS_values.txt > temp.csv
	cat temp.csv >> OMEGA_per_branch_genetree/Opsins.omega.csv

	grep "LB\":" Reconciled_trees_95/$clade_name/Omega_per_branch/$clade_name.fish.cds.trimmed.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
	grep "MLE\":" Reconciled_trees_95/$clade_name/Omega_per_branch/$clade_name.fish.cds.trimmed.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
	grep "UB\":" Reconciled_trees_95/$clade_name/Omega_per_branch/$clade_name.fish.cds.trimmed.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
	grep "\"dN\"" Reconciled_trees_95/$clade_name/Omega_per_branch/$clade_name.fish.cds.trimmed.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
	grep "\"dS\"" Reconciled_trees_95/$clade_name/Omega_per_branch/$clade_name.fish.cds.trimmed.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
	grep -B2 "LB\":" Reconciled_trees_95/$clade_name/Omega_per_branch/$clade_name.fish.cds.trimmed.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels
	sed -i "s/^/$clade_name,/g" curr_labels


	paste -d "," curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_labels curr_dN_values.txt curr_dS_values.txt > temp.csv
	cat temp.csv >> OMEGA_per_branch_speciestree/Opsins.omega.csv

done






for clade in Clades/* ; do 

	clade_name=`echo "$clade" | sed 's/Clades\///g'`

	sed "s/^/$clade_name,/g" Reconciled_trees_95/$clade_name/$clade_name.renaming_file.csv >> OMEGA_per_branch_speciestree/Corresponding_labels.csv 

done






#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
###################### Accessory scripts ##############################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################


====== run_fitmg94.sh ===========================================================================================

#!/bin/bash


#SBATCH --job-name=fitMG   # Job name

eval "$(conda shell.bash hook)"
conda activate HyPhy_Env


clade_name=$1


hyphy FitMG94.bf --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile --code Universal --type local ENV=TOLERATE_NUMERICAL_ERRORS=1


======================================================================================================================================================================================




============================================ iqtree_phylo.sh 


#!/bin/bash


#SBATCH --job-name=Iqtree_opsin   # Job name

module load IQ-TREE/2.0-rc1-foss-2018b

gene_alignment=$1 

iqtree -s $gene_alignment --seqtype AA -nt 20 -bb 1000 


======================================================================================================================================================================================




============================================ iqtree_all_phylo.sh 


#!/bin/bash


#SBATCH --job-name=Iqtree_opsin   # Job name

module load IQ-TREE/2.0-rc1-foss-2018b

gene_alignment=$1 

iqtree -s $gene_alignment --seqtype AA -m JTT+F+G4 -nt 40 -bb 1000 

======================================================================================================================================================================================



============================================ root_opsin_phylo.R

#load packages

library("ape")
library("dplyr")
library("phytools")
library("ggtree")

#load the tree
mytree <- read.tree("Complete_opsins.aln.trimal.treefile")
#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

write.tree(mytree_label, "Complete_opsins.aln.trimal.treefile.nodelabel")
#Import the ID of one exorh and one opn4x gene
tip1="Plectropomus_leopardus---opn7b-CM018415.2-28786147-28789316---5_exons"
tip2="Hoplias_malabaricus---opn4x-CM056246.1-9517823-9522094---7_exons"


#Root the tree at the common ancestor of opn4x and exorh
MRCA_opsins <- findMRCA(mytree_label, tips=c(tip1, tip2), type="node") #extract the node name
mytree_rooted <- root(mytree_label, node=MRCA_opsins, resolve.root= TRUE) #root


#write in file
write.tree(mytree_rooted, file="Complete_opsins.aln.trimal.treefile.nodelabel.rooted")


======================================================================================================================================================================================



============================================ Remove_outgroup.R


library(ape)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
outgroup <- scan(args[2], what="character")

mytree_pruned <- drop.tip(mytree, outgroup)

write.tree(mytree_pruned, file=args[3])


======================================================================================================================================================================================


============================================ Collapse_tree_90.R


library(ape)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
Badnodes <- which(as.numeric(mytree$node.label) < 90) + length(mytree$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){ Badnodes_indexes <- c(Badnodes_indexes, which(mytree$edge[,2] == node)) }

mytree$edge.length[Badnodes_indexes] <- 0 
tree_multi <- di2multi(mytree) 
write.tree(tree_multi, file = args[2])


======================================================================================================================================================================================



============================================ Collapse_tree_95.R


library(ape)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
Badnodes <- which(as.numeric(mytree$node.label) < 95) + length(mytree$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){ Badnodes_indexes <- c(Badnodes_indexes, which(mytree$edge[,2] == node)) }

mytree$edge.length[Badnodes_indexes] <- 0 
tree_multi <- di2multi(mytree) 
write.tree(tree_multi, file = args[2])



======================================================================================================================================================================================




====== run_hyphy_seq.sh ===========================================================================================

#!/bin/bash


#SBATCH --job-name=Multip_hphy   # Job name


module load HyPhy/2.5.31-foss-2018b 


clade_name=$1


hyphy slac --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile --code Universal
hyphy fel --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile --code Universal
hyphy absrel --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile --code Universal


======================================================================================================================================================================================

====== run_absrel_seq.sh ===========================================================================================

#!/bin/bash


#SBATCH --job-name=absrel   # Job name


module load HyPhy/2.5.31-foss-2018b 


clade_name=$1

hyphy absrel --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile --code Universal



======================================================================================================================================================================================


====== run_bgm_seq.sh ===========================================================================================


#!/bin/bash


#SBATCH --job-name=bgm   # Job name


module load HyPhy/2.5.31-foss-2018b 

clade_name=$1

hyphy bgm --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile --code Universal


======================================================================================================================================================================================




============================================ run_busted.sh 


#!/bin/bash

module load HyPhy/2.5.31-foss-2018b 

clade_name=$1 

hyphy busted --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.myphylo_noloss.nwk --code Universal --branches TestBranch



======================================================================================================================================================================================



============================================ run_FUBAR_sptree.sh 


#!/bin/bash

module load HyPhy/2.5.31-foss-2018b 

clade_name=$1 

hyphy fubar --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.myphylo_noloss.nwk


======================================================================================================================================================================================



============================================ run_FUBAR.sh 


#!/bin/bash

module load HyPhy/2.5.31-foss-2018b 

clade_name=$1 

hyphy fubar --alignment $clade_name.fish.cds.trimmed.aln --tree $clade_name.prot.trimmed.aln.nooutgroup.treefile



======================================================================================================================================================================================






====== prepare_codeml_dup.R ===========================================================================================



library("caper")
library("ape")
library(adephylo)
library("phytools")
library("aphylo")
library("treeio")
library(dplyr)
library(data.table)
library(tidyverse)
library("phangorn")

args = commandArgs(trailingOnly=TRUE)

mytree <- 
  read.nhx(args[1])

#rename Nodes

mytree@phylo <- makeNodeLabel(mytree@phylo, method = "number", prefix = "Node")

species_test <- "test"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(mytree@phylo, misc_table, by = 'node') %>%
  dplyr::select(node, label)


#Extract branches resulting from a duplication

dup_nodes <- 
  mytree@data %>%
  filter(D == "Y") %>%
  pull(node)

daughter_dups_all <- c()
for(mydup in dup_nodes){
  
  daughter_dups <- 
    Children(mytree@phylo,mydup)
  
  daughter_dups_all <- 
    c(daughter_dups_all,
      daughter_dups)
  
}

daughter_dups_all_label <- 
  node_label_corresp %>%
  filter(node %in% daughter_dups_all) %>%
  pull(label)


#remove LOST branches

myphylo <- mytree@phylo
LOST_branches <- 
  node_label_corresp[grep("LOST", node_label_corresp$label), ] %>%
  pull(label)
myphylo_noloss <- 
  ape::drop.tip(myphylo, LOST_branches)

#remove branch length
myphylo_noloss$edge.length <- NULL


#write the tree and alignment for PAML

#=> write(daughter_dups_all_label, "~/Non_visual_opsins_Project/Phylogenies_Opsins_v2/Reconciled_trees_95/dup_labels.txt")
write(daughter_dups_all_label, args[2])


#add a "#1" for PAML to nodes

all_tips_labels <- myphylo_noloss$tip.label
#=> write(all_tips_labels, "~/Non_visual_opsins_Project/Phylogenies_Opsins_v2/Reconciled_trees_95/all_tips_labels")
write(all_tips_labels, args[3])


for(curr_tip in all_tips_labels){
  
  if(curr_tip %in% daughter_dups_all_label){
    
    curr_tip_modif <- paste(curr_tip, " #1", sep="")
    all_tips_labels[all_tips_labels == curr_tip] <- curr_tip_modif
    
  }
}
myphylo_noloss$tip.label <- all_tips_labels 



all_node_labels <- myphylo_noloss$node.label
#=> write(all_node_labels, "~/Non_visual_opsins_Project/Phylogenies_Opsins_v2/Reconciled_trees_95/all_node_labels")
write(all_node_labels, args[4])


for(curr_node in all_node_labels){
  
  if(curr_node %in% daughter_dups_all_label){
    
    curr_node_modif <- paste(curr_node, " #1", sep="")
    all_node_labels[all_node_labels == curr_node] <- curr_node_modif
    
  }
}
myphylo_noloss$node.label <- all_node_labels 

#=> write.tree(myphylo_noloss, "~/Non_visual_opsins_Project/Phylogenies_Opsins_v2/Reconciled_trees_95/myphylo_noloss.nwk")
write.tree(myphylo_noloss, args[5])


======================================================================================================================================================================================


