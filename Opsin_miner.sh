#!/bin/bash


#SBATCH --job-name=Opsin_mining   # Job name

#eval "$(conda shell.bash hook)"
#conda activate olfactory

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module load StdEnv/2020
module load gcc/9.3.0
module load r/4.2.1
module load blast+/2.12.0
module load emboss/6.6.0
module load samtools/1.15.1
module load mafft/7.471
module load iq-tree/2.0.7
module load python/3.9.6
module load fastx-toolkit/0.0.14

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#Initialize variables 

genome=$1 #genome as fasta file
opsin_database=$2 #database of opsins as fasta (proteins)
blast_database=$3 #blast database (uniprot db, proteins)
scripts_folder_location=$4 ; scripts_location=`echo "$scripts_folder_location" | sed 's/\/$//'` #scripts database
maximum_intron_length=$5 #maximum intron length : will define the extension length
number_of_thread=$6 #number of threads to use
maximum_intron_length_half=$((maximum_intron_length / 2)) 



#Makeblastdb so we can blast genes against the genome

if test -f "$genome.ndb" ; then echo "Genome blast database already exist" ; else makeblastdb -in $genome -dbtype nucl ; fi 
if test -f "$genome.fai" ; then echo "Genome fai file already exist" ; else samtools faidx $genome ; fi 
echo $genome
echo $opsin_database
#Perform tblastn using known opsin genes against the genome with an evalue of 1e-05

tblastn -query $opsin_database -db $genome -evalue 1e-05 -outfmt 6 -out Opsin_vs_Genome.blastn -num_threads $number_of_thread

#Lets launch a Rscript that will merge all blast hits 

Rscript $scripts_location/Merge_blast_hits.R 
xargs samtools faidx $genome < Blast_nonoverlapping.tsv > Blast_nonoverlapping.fasta

#Retain only besthit that best match to an opsin gene
blastx -query Blast_nonoverlapping.fasta -db Database/GPCR_plus_Olfactory_plus_Taste_receptors_vertebrates_reformat.prot -max_target_seqs 1 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_thread
grep "Opsin-" blastx_blast_regions.tsv | cut -f1 | sort | uniq > Opsins_best_hits.txt
for i in `cat Opsins_best_hits.txt` ; do grep "$i" Blast_nonoverlapping.tsv >> Opsins_Regions.tsv ; done


#Extend all best hits by Xbp (X=intron size) upstream and downstream . Result file : Potential_Opsins_regions.tsv
Rscript $scripts_location/Merge_blast_hits_extend.R $maximum_intron_length



#Split the opsin database and launch exonerate with these sequences against potential opsin regions (max intron length : 30000bp)

mkdir Splitted_db
$scripts_location/exonerate-2.2.0-x86_64/bin/fastasplit -f Database/all_opsins.prot -c 30 --output Splitted_db



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Search for Opsins genes in a loop  #############################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


mkdir COFFRE_PREDICTIONS/


#re-initialize files

[ -e Potential_multiple_exon_CDS.fa ] && rm Potential_multiple_exon_CDS.fa 
[ -e Pseudogenes_multiple_exon.fa ] && rm Pseudogenes_multiple_exon.fa 
[ -e No_opsins_genes_coordinates.txt ] && rm No_opsins_genes_coordinates.txt
[ -e Frameshift_less_Pseudogenes.fa ] && rm Frameshift_less_Pseudogenes.fa

#Start the loop to search for Opsins genes

current_nb_sequences=1
previous_iteration_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
number_regions_blast=`grep "[0-9]" Potential_Opsins_regions.tsv | wc -l`


while [ "$current_nb_sequences" -gt "$previous_iteration_nb_sequences" ] && [ "$number_regions_blast" -gt "0" ] ; do


	previous_iteration_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
	
	#Extract identified regions in a fasta file

	xargs samtools faidx $genome < Potential_Opsins_regions.tsv > Potential_Opsins_regions.fa

	blastx -query Potential_Opsins_regions.fa -db Database/all_opsins.prot -max_target_seqs 5 -outfmt '6 qseqid sseqid' -out blastx_per_scaffold.tsv -num_threads $number_of_thread
	grep ">" Potential_Opsins_regions.fa | sed 's/>//g' > scaffold_id.txt
	
	rm -r Gene_to_exo_per_scaff ; mkdir Gene_to_exo_per_scaff
	for i in `cat scaffold_id.txt` ; do grep "$i" blastx_per_scaffold.tsv | cut -f2 | sed 's/-.*//g' | sort | uniq > Gene_to_exo_per_scaff/$i.target ; done
	for file in Gene_to_exo_per_scaff/* ; do for gene in `cat $file` ; do cat Database/Final_dataset/$gene.prot >> $file.prot ; done ; done
	
	
	for file in Gene_to_exo_per_scaff/*.target ; do scaffold_name=`echo "$file" | sed 's/.target//g' | sed 's/Gene_to_exo_per_scaff\///g'` ; samtools faidx $genome $scaffold_name > Gene_to_exo_per_scaff/$scaffold_name.fasta  ; done
	
	for file in Gene_to_exo_per_scaff/* ; do file_name_pars=`echo "$file" | sed 's/:/-/g'` ; mv $file $file_name_pars ; done


	mkdir Exonerate_raw_results_folder
	
	for i in Gene_to_exo_per_scaff/*.fasta ; do
		file_name=`echo $i | sed 's/Gene_to_exo_per_scaff\///g' | sed 's/.fasta//g'`
		prot_to_exo=`echo $file_name.target.prot`
		sbatch --account=def-mshafer --job-name=ops_vs_scaff -W -c 2 --qos=6hours --mem=4G --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.exo.rslt ; sleep 10" &
		sbatch --account=def-mshafer --job-name=ops_vs_scaff -W -c 2 --qos=6hours --mem=4G --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.noexhaustive.exo.rslt ; sleep 10" &
	done


	echo "Exonerate running -- Wait"
	
	wait
	
	echo "Exonerate research done"

	
	#Merge exonerate results
	
	cat Exonerate_raw_results_folder/*.exo.rslt > Exonerate_results.txt
	
	
	#extract vulgar lines
	
	grep "vulgar" Exonerate_results.txt > vulgar_lines.txt 
	
	
	#extract interesting columns of vulgar lines
	#query, query_start, query_end, scaffold, scaffold_start, scaffold_end, strand, exonerate_score
	
	sed 's/vulgar: //g' vulgar_lines.txt | cut -f1,2,3,5,6,7,8,9 -d " " > vulgar_lines_parsed.txt
	
	
	#count the number of introns using vulgar lines
	
	
	IFS=$'\n'
	
	awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/ I /,"") "\t" NR}' vulgar_lines.txt > number_introns_per_line.txt
	grep -v "count" number_introns_per_line.txt | cut -f1  > intron_numbers.txt 
	
	#add the intron number to vulgar lines 
	
	paste -d " " vulgar_lines_parsed.txt intron_numbers.txt > vulgar_lines_intron_numbers.txt
	
	
	##Add informations about the best blastp results of each exonerate predicted genes
	
	
	sed -n '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/p' Exonerate_results.txt | sed 's/C4 Alignment:.*//g' | sed 's/Hostname:.*//g' | sed 's/Command line:.*//g' | sed 's/^--$//g' | sed 's/-- completed exonerate analysis.*//g' | sed 's/# --- END OF GFF DUMP ---//g' | sed 's/^#$/>seq_to_rename/g' > List_exonerate_cds.fasta #extract all predicted genes sequences
	transeq List_exonerate_cds.fasta List_exonerate_cds.prot #translate sequences
	sed 's/ /_/g' vulgar_lines_intron_numbers.txt > sequences_names.txt #extract exonerate vulgar line to rename sequences
	awk '/^>/ { printf("%s_%s\n",$0,i++);next;} { print $0;}' List_exonerate_cds.prot > List_exonerate_cds_renamed.prot #first round of rename
	awk '/^>/ { printf("%s_%s\n",$0,i++);next;} { print $0;}' List_exonerate_cds.fasta > List_exonerate_cds_renamed.fasta
	grep ">" List_exonerate_cds_renamed.prot | sed 's/>//g' > old_names.txt #extract names
	paste -d "\t" old_names.txt sequences_names.txt > renaming_file #crate a file for rename_fasta.pl
	perl $scripts_location/rename_fasta.pl renaming_file List_exonerate_cds_renamed.prot > List_exonerate_cds.prot #completely rename sequences with the exonerate vulgar line
	perl $scripts_location/rename_fasta.pl renaming_file List_exonerate_cds_renamed.fasta > List_exonerate_cds.fasta

	#Perform the blastp
	blastx -query List_exonerate_cds.fasta -db Database/all_opsins.prot -outfmt '6 qseqid sseqid evalue' -out all_blastp.txt -max_target_seqs 1 -num_threads $number_of_thread
	
	#Extract the information
	[ -e all_blastp_parsed.txt ] && rm all_blastp_parsed.txt
	for i in `cat sequences_names.txt` ; do if grep -q "$i" all_blastp.txt ; then grep -m1 "$i" all_blastp.txt | cut -f2,3 >> all_blastp_parsed.txt ; else echo "NoQuery	99999" >> all_blastp_parsed.txt ; fi ; done 
	sed -i 's/	/ /g' all_blastp_parsed.txt
	paste -d " " vulgar_lines_intron_numbers.txt all_blastp_parsed.txt > vulgar_lines_intron_numbers_blastrslt.txt
	
	
	
	#Parse exonerate results. Find the best exonerate results that are most likely complete genes or pseudoogenes, and not overlapping
	
	Rscript $scripts_location/Parse_exonerate_results.R #result file : Parsed_exonerate_gene_regions.tsv
	

	nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`
	if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then


		IFS=$'\n'
		
		for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
			
			query=`echo "$line" | cut -f7`
			scaffold=`echo "$line" | cut -f1`
			scaff_start=`echo "$line" | cut -f2`
			scaff_end=`echo "$line" | cut -f3`
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv
		
		
		#Let's now predict genes on these regions !
		
		IFS=$'\n'
		
		
		mkdir Genes_predictions
		

		#for line in `cat Correct_coordinates_for_exonerate.tsv` ; do scaffold_s_e=`echo "$line" | cut -f1` ; best_query=`echo "$line" | cut -f2` ; scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'` ; samtools faidx $genome $scaffold_s_e > scaffold.fa ; sed -i 's/:/-/g' scaffold.fa ; samtools faidx Database/all_opsins.prot $best_query > query.prot ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate ; cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta ; if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos ; sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence ; cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate ; fi ; done
		
		for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
			
			scaffold_s_e=`echo "$line" | cut -f1`
			best_query=`echo "$line" | cut -f2`
			scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
		
		
			samtools faidx $genome $scaffold_s_e > scaffold.fa
			sed -i 's/:/-/g' scaffold.fa
		
			samtools faidx Database/all_opsins.prot $best_query > query.prot
		
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
		
			#if exonerate with E=TRUE failed, then launch a normal exonerate

			exonerate_output_linenb=`wc -l < Genes_predictions/$scaffold_s_e_n.exonerate`
			if [ "$exonerate_output_linenb" -lt "10" ] ; then
				$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron 50000 --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			fi


			#extract only the best result if there are two with the same score
		
			if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
				sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
				sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
				cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
			fi
		
		
		
		done
		
		
		### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
		
		
		#Result folder
		mkdir Filtered_predictions	
		
		for file in Genes_predictions/*.exonerate ; do
		
			#extract some infos from file name
			file_name=`echo "$file" | sed 's/.*\///g'`
			file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
			fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
			initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
		
			#Test if the predicted gene is a Opsins gene or not
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot
			blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_thread
			
		
			#Lets continue only if the best match is an OR
			#if grep -q -i "olfactory\|odorant" blastp_result ; then 
			if grep -q -i "RGR_\|RPE-retinal\|opsin" blastp_result ; then
		
				#Define the scaffold  name
				scaffold=`echo "$file" | sed 's/.*\///g' | sed 's/-.*//g'`
				#Define the strand on which the predicted gene is
				strand=`grep "	similarity	" $file | cut -f7`
				#Define the first position of the query on the target sequence
				first_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
				#Define the last position of the query on the target sequence
				second_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
				
				
				#Lets extract CDS if the gene is on the negative strand
				if [ $strand == "-" ] ; then 
				
				#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
		
		
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
		
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
					query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
					perc80_query_length=$((query_length*80/100*3))
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3
					if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
				
		
		
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


						extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
		
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl
						samtools faidx Filtered_predictions/$file_name_reduced.ORFP
						gene_length=`cut -f2 Filtered_predictions/$file_name_reduced.ORFP.fai | head -1`
						rm Filtered_predictions/$file_name_reduced.ORFP.fai
						tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 > current_verification_exons.tsv
						number_blast_hit=`awk '{ if ($4 >= 60) { print } }' current_verification_exons.tsv | wc -l`
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
						length_align=`awk '{ if ($4 >= 30) { print } }' current_verification_exons.tsv | cut -f7,8 | awk '{ $3 = $2 - $1 } 1' | cut -f3 -d " " | awk -F',' '{sum+=$1;} END{print sum;}'`
						half_gene_length=$((gene_length / 2))
						gene_length_plusfiftyperc=$((gene_length + half_gene_length))
						if [ "$length_align" -gt "$gene_length_plusfiftyperc" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 

		
		
					#If not ORF found, then determinate the gene state
				
					elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
				
						stop_codon_state="FALSE"
						edge_state="FALSE"
						frameshift_state="FALSE"
				
						##Stop codon checking
				
						#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
						transeq predicted_cds.fa predicted_cds.prot
				
						#Estimate the interval on which we wil search stop codons. 
						query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
						query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
						query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
				
						#Lets see if we find stop codon before the five_percent_position
						stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
				
						if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
				
				
						##Frameshift checking
				
						#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
						#We remove border exons if there are less than 60nt in length. Run as iteration.
				
						grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
						awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
				
						#Check for the presence of frameshift
						frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
						if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
						if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
				
				
				
						##Edge checking
				
						#Check if the gene is at a conting border
						#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
						gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
						gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
						true_start_coord=$((extracted_scaffold_start + gene_start_coord))
						true_end_coord=$((extracted_scaffold_start + gene_end_coord))
				
						#First check if these coordinates are near the end of scaffolds (<5000 bp)
						if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
						scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
						diff_lengths=$((scaffold_length - true_end_coord))
						if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
						
						#Now check if there are consecutive N near the gene that could indicate conting end
						extanded_start_coord=$((true_start_coord - 200))
						extanded_end_coord=$((true_end_coord + 200))
				
						#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
						consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
						if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
						if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				
						##Extract the sequence
				
						[ -e Current_exon_rev.txt ] && rm Current_exon_rev.txt
				
						#Extract the corresponding sequence
						for line in `cat Correct_exons.txt` ; do
							start_pos=`echo "$line" | cut -f1`
							end_pos=`echo "$line" | cut -f2`
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
							revseq Current_exon.fa Current_exon_rev.fa
				
							
				
							#add the reversed sequence to a text file
							grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
				
						done
				
						#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
						exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
				
						header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
						sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
						sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
						cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
						sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
	
	
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl
						samtools faidx predicted_cds.prot
						gene_length=`cut -f2 predicted_cds.prot.fai | head -1`
						rm predicted_cds.prot.fai

						tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 > current_verification_exons.tsv
						number_blast_hit=`awk '{ if ($4 >= 60) { print } }' current_verification_exons.tsv | wc -l`
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
						length_align=`awk '{ if ($4 >= 30) { print } }' current_verification_exons.tsv | cut -f7,8 | awk '{ $3 = $2 - $1 } 1' | cut -f3 -d " " | awk -F',' '{sum+=$1;} END{print sum;}'`
						half_gene_length=$((gene_length / 2))
						gene_length_plusfiftyperc=$((gene_length + half_gene_length))
						if [ "$length_align" -gt "$gene_length_plusfiftyperc" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$length_align" -gt "$gene_length_plusfiftyperc" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 

				
					fi
		
		
		
				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
					query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
					perc80_query_length=$((query_length*80/100*3))
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3 -reverse FALSE
				
				
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


						extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
					
		
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl
						samtools faidx Filtered_predictions/$file_name_reduced.ORFP
						gene_length=`cut -f2 Filtered_predictions/$file_name_reduced.ORFP.fai | head -1`
						rm Filtered_predictions/$file_name_reduced.ORFP.fai
						tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 > current_verification_exons.tsv
						number_blast_hit=`awk '{ if ($4 >= 60) { print } }' current_verification_exons.tsv | wc -l`
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
						length_align=`awk '{ if ($4 >= 30) { print } }' current_verification_exons.tsv | cut -f7,8 | awk '{ $3 = $2 - $1 } 1' | cut -f3 -d " " | awk -F',' '{sum+=$1;} END{print sum;}'`
						half_gene_length=$((gene_length / 2))
						gene_length_plusfiftyperc=$((gene_length + half_gene_length))
						if [ "$length_align" -gt "$gene_length_plusfiftyperc" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 

		
					#If not ORF found, then determinate the gene state
				
					elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
				
						stop_codon_state="FALSE"
						edge_state="FALSE"
						frameshift_state="FALSE"
				
					##Stop codon checking
				
						#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
						transeq predicted_cds.fa predicted_cds.prot
				
						#Estimate the interval on which we wil search stop codons. 
						query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
						query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
						query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
				
						#Lets see if we find stop codon before the five_percent_position
						stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
				
						if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
				
				
						##Frameshift checking
				
						#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
						#We remove border exons if there are less than 60nt in length. Run as iteration.
				
						grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
						awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
				
						#Check for the presence of frameshift
						frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
						if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
						if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
				
				
				
						##Edge checking
				
						#Check if the gene is at a conting border
						#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
						gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
						gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
						true_start_coord=$((extracted_scaffold_start + gene_start_coord))
						true_end_coord=$((extracted_scaffold_start + gene_end_coord))
				
						#First check if these coordinates are near the end of scaffolds (<5000 bp)
						if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
						scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
						diff_lengths=$((scaffold_length - true_end_coord))
						if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
						
						#Now check if there are consecutive N near the gene that could indicate conting end
						extanded_start_coord=$((true_start_coord - 200))
						extanded_end_coord=$((true_end_coord + 200))
				
						#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
						consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
						if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
						if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				
						##Extract the sequence
				
						[ -e Current_exon.txt ] && rm Current_exon.txt
						#Extract the corresponding sequence
						for line in `cat Correct_exons.txt` ; do
							start_pos=`echo "$line" | cut -f1`
							end_pos=`echo "$line" | cut -f2`
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
				
							#add the reversed sequence to a text file
							grep -v ">" Current_exon.fa >> Current_exon.txt
				
						done
				
				
						#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
						exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
				
						header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
						sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
						sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
						cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
						sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl
						samtools faidx predicted_cds.prot
						gene_length=`cut -f2 predicted_cds.prot.fai | head -1`
						rm predicted_cds.prot.fai

						tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 > current_verification_exons.tsv
						number_blast_hit=`awk '{ if ($4 >= 60) { print } }' current_verification_exons.tsv | wc -l`
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
						length_align=`awk '{ if ($4 >= 30) { print } }' current_verification_exons.tsv | cut -f7,8 | awk '{ $3 = $2 - $1 } 1' | cut -f3 -d " " | awk -F',' '{sum+=$1;} END{print sum;}'`
						half_gene_length=$((gene_length / 2))
						gene_length_plusfiftyperc=$((gene_length + half_gene_length))
						if [ "$length_align" -gt "$gene_length_plusfiftyperc" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$length_align" -gt "$gene_length_plusfiftyperc" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 

				
	
	
					fi
				fi
		
			else echo "$initial_header" >> No_opsins_genes_coordinates.txt
		
			fi
		
		done
	fi
	
	
	#Now that we have filtered all our results, we can concatenate the results
	
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa
 

 	#remove redundant sequences names
 	awk '/^>/{f=!d[$1];d[$1]=1}f' Potential_multiple_exon_CDS.fa > non_redundant.fa ; mv non_redundant.fa Potential_multiple_exon_CDS.fa


	#Extract coordinates of found genes
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	sed 's/-/	/g' No_opsins_genes_coordinates.txt >> Coordinates_already_examined.tsv
	if [ `wc -l < Coordinates_already_examined.tsv` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi
	
	current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
	
	#re-process blast result to find potential Opsins regions exclusing already found genes
	
	Rscript $scripts_location/Merge_blast_hits_extend_filterAlreadyFound.R $maximum_intron_length


	cp Filtered_predictions/* COFFRE_PREDICTIONS/
	cp Genes_predictions/* COFFRE_PREDICTIONS/

	rm -r Filtered_predictions/
	rm -r Genes_predictions/
	rm -r Exonerate_raw_results_folder/
	rm Parsed_exonerate_gene_regions.tsv
	

	if test -f "Potential_Opsins_regions.tsv" ; then number_regions_blast=`grep "[0-9]" Potential_Opsins_regions.tsv | wc -l` ; else number_regions_blast=0 ; fi

	
done


rm -r Filtered_predictions/
rm -r Genes_predictions/
rm Parsed_exonerate_gene_regions.tsv



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Search for remaining Opsins genes with size a bit below  #######################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


Rscript $scripts_location/Parse_exonerate_results_second.R

nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`


if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then

	IFS=$'\n'
	
	for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
		
		query=`echo "$line" | cut -f7`
		scaffold=`echo "$line" | cut -f1`
		scaff_start=`echo "$line" | cut -f2`
		scaff_end=`echo "$line" | cut -f3`
	
		echo "$scaffold:$scaff_start-$scaff_end	$query"
	
	done > Correct_coordinates_for_exonerate.tsv
	
	
	#Let's now predict genes on these regions !
	
	IFS=$'\n'
	
	
	mkdir Genes_predictions
	
	
	for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
		
		scaffold_s_e=`echo "$line" | cut -f1`
		best_query=`echo "$line" | cut -f2`
		scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
	
	
		samtools faidx $genome $scaffold_s_e > scaffold.fa
		sed -i 's/:/-/g' scaffold.fa
	
		samtools faidx Database/all_opsins.prot $best_query > query.prot
	
		$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
		cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
	

		#if exonerate with E=TRUE failed, then launch a normal exonerate

		exonerate_output_linenb=`wc -l < Genes_predictions/$scaffold_s_e_n.exonerate`
		if [ "$exonerate_output_linenb" -lt "10" ] ; then
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron 50000 --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
		fi

	
		#extract only the best result if there are two with the same score
	
		if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
			sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
			sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
			cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
		fi
	
	
	
	done
	
	
	### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
	
	
	#Result folder
	mkdir Filtered_predictions	
	
	for file in Genes_predictions/*.exonerate ; do
	
		#extract some infos from file name
		file_name=`echo "$file" | sed 's/.*\///g'`
		file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
		fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
		initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
	
		#Test if the predicted gene is a Opsins gene or not
		awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
		transeq predicted_cds.fa predicted_cds.prot
		blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
		
	
		#Lets continue only if the best match is an OR
		#if grep -q -i "olfactory\|odorant" blastp_result ; then 
		if grep -q -i "RGR_\|RPE-retinal\|opsin" blastp_result ; then
	
			#Define the scaffold  name
			scaffold=`echo "$file" | sed 's/.*\///g' | sed 's/-.*//g'`
			#Define the strand on which the predicted gene is
			strand=`grep "	similarity	" $file | cut -f7`
			#Define the first position of the query on the target sequence
			first_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			#Define the last position of the query on the target sequence
			second_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			
			
			#Lets extract CDS if the gene is on the negative strand
			if [ $strand == "-" ] ; then 
			
			#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
	
	
				#If strand is minus, then the first position is:
				target_end=$((first_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((first_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
			
				#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
	
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
				query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
				perc80_query_length=$((query_length*80/100*3))
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3
				if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
			
	
	
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
			
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 60) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
	
	
	
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
					##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<5000 bp)
					if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon_rev.txt ] && rm Current_exon_rev.txt
			
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
						revseq Current_exon.fa Current_exon_rev.fa
			
						
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
			
					done
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
	
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 60) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
	
	
			
				fi
	
	
	
			#Lets make the same steps with slight modifications for the + strand
			elif [ $strand == "+" ] ; then 
			
				#If strand is minus, then the first position is:
				target_end=$((second_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((second_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
			
				#Extract the CDS sequence predicted by exonerate and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
				query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
				perc80_query_length=$((query_length*80/100*3))
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3 -reverse FALSE
			
			
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 60) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
	
	
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
				##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<5000 bp)
					if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon.txt ] && rm Current_exon.txt
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon.fa >> Current_exon.txt
			
					done
			
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 60) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
	
	
				fi
			fi
	
		else echo "$initial_header" >> No_opsins_genes_coordinates.txt
	
		fi
	
	done

fi

#Now that we have filtered all our results, we can concatenate the results

for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa

#Extract coordinates of found genes
grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
sed 's/-/	/g' No_opsins_genes_coordinates.txt >> Coordinates_already_examined.tsv
if [ `wc -l < Coordinates_already_examined.tsv` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi

current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`

#re-process blast result to find potential Opsins regions exclusing already found genes


cp Filtered_predictions/* COFFRE_PREDICTIONS/
cp Genes_predictions/* COFFRE_PREDICTIONS/

rm -r Filtered_predictions/
rm -r Genes_predictions/
rm Parsed_exonerate_gene_regions.tsv


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Search for Opsins pseudogenes  #################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

Rscript $scripts_location/Parse_exonerate_results_third.R



nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`


if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then

	IFS=$'\n'
	
	for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
		
		query=`echo "$line" | cut -f7`
		scaffold=`echo "$line" | cut -f1`
		scaff_start=`echo "$line" | cut -f2`
		scaff_end=`echo "$line" | cut -f3`
	
		echo "$scaffold:$scaff_start-$scaff_end	$query"
	
	done > Correct_coordinates_for_exonerate.tsv
	
	
	#Let's now predict genes on these regions !
	
	IFS=$'\n'
	
	
	mkdir Genes_predictions
	
	
	for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
		
		scaffold_s_e=`echo "$line" | cut -f1`
		best_query=`echo "$line" | cut -f2`
		scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
	
	
		samtools faidx $genome $scaffold_s_e > scaffold.fa
		sed -i 's/:/-/g' scaffold.fa
	
		samtools faidx Database/all_opsins.prot $best_query > query.prot
	
		$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
		cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
	
		#if exonerate with E=TRUE failed, then launch a normal exonerate

		exonerate_output_linenb=`wc -l < Genes_predictions/$scaffold_s_e_n.exonerate`
		if [ "$exonerate_output_linenb" -lt "10" ] ; then
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron 50000 --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
		fi

	
		#extract only the best result if there are two with the same score
	
		if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
			sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
			sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
			cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
		fi
	
	
	
	done
	
	
	### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
	
	
	#Result folder
	mkdir Filtered_predictions	
	
	for file in Genes_predictions/*.exonerate ; do
	
		#extract some infos from file name
		file_name=`echo "$file" | sed 's/.*\///g'`
		file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
		fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
		initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
	
		#Test if the predicted gene is a Opsins gene or not
		awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
		transeq predicted_cds.fa predicted_cds.prot
		blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
		
	
		#Lets continue only if the best match is an OR
		#if grep -q -i "olfactory\|odorant" blastp_result ; then 
		if grep -q -i "RGR_\|RPE-retinal\|opsin" blastp_result ; then
	
			#Define the scaffold  name
			scaffold=`echo "$file" | sed 's/.*\///g' | sed 's/-.*//g'`
			#Define the strand on which the predicted gene is
			strand=`grep "	similarity	" $file | cut -f7`
			#Define the first position of the query on the target sequence
			first_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			#Define the last position of the query on the target sequence
			second_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			
			
			#Lets extract CDS if the gene is on the negative strand
			if [ $strand == "-" ] ; then 
			
			#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
	
	
				#If strand is minus, then the first position is:
				target_end=$((first_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((first_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
			
				#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
	
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
				query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
				perc80_query_length=$((query_length*80/100*3))
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3
				if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
			
	
	
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
			
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
					##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<5000 bp)
					if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon_rev.txt ] && rm Current_exon_rev.txt
			
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
						revseq Current_exon.fa Current_exon_rev.fa
			
						
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
			
					done
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
	
			
				fi
	
	
	
			#Lets make the same steps with slight modifications for the + strand
			elif [ $strand == "+" ] ; then 
			
				#If strand is minus, then the first position is:
				target_end=$((second_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((second_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
			
				#Extract the CDS sequence predicted by exonerate and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
				query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
				perc80_query_length=$((query_length*80/100*3))
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3 -reverse FALSE
			
			
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
				##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<5000 bp)
					if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon.txt ] && rm Current_exon.txt
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon.fa >> Current_exon.txt
			
					done
			
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
				fi
			fi
	
		else echo "$initial_header" >> No_opsins_genes_coordinates.txt
	
		fi
	
	done
fi


#Now that we have filtered all our results, we can concatenate the results

for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa

#Extract coordinates of found genes
grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
sed 's/-/	/g' No_opsins_genes_coordinates.txt >> Coordinates_already_examined.tsv
if [ `wc -l < Coordinates_already_examined.tsv` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi

current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`

#re-process blast result to find potential Opsins regions exclusing already found genes


cp Filtered_predictions/* COFFRE_PREDICTIONS/
cp Genes_predictions/* COFFRE_PREDICTIONS/


rm -r Filtered_predictions/
rm -r Genes_predictions/
rm Parsed_exonerate_gene_regions.tsv



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Final round to find Opsins pseudogenes with length below previous iteration  ###################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


Rscript $scripts_location/Parse_exonerate_results_final.R

nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`


if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then	

	IFS=$'\n'
	
	for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
		
		query=`echo "$line" | cut -f7`
		scaffold=`echo "$line" | cut -f1`
		scaff_start=`echo "$line" | cut -f2`
		scaff_end=`echo "$line" | cut -f3`
	
		echo "$scaffold:$scaff_start-$scaff_end	$query"
	
	done > Correct_coordinates_for_exonerate.tsv
	
	
	#Let's now predict genes on these regions !
	
	IFS=$'\n'
	
	
	mkdir Genes_predictions
	
	
	for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
		
		scaffold_s_e=`echo "$line" | cut -f1`
		best_query=`echo "$line" | cut -f2`
		scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
	
	
		samtools faidx $genome $scaffold_s_e > scaffold.fa
		sed -i 's/:/-/g' scaffold.fa
	
		samtools faidx Database/all_opsins.prot $best_query > query.prot
	
		$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
		cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
	

		#if exonerate with E=TRUE failed, then launch a normal exonerate

		exonerate_output_linenb=`wc -l < Genes_predictions/$scaffold_s_e_n.exonerate`
		if [ "$exonerate_output_linenb" -lt "10" ] ; then
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron 50000 --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
		fi


		#extract only the best result if there are two with the same score
	
		if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
			sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
			sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
			cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
		fi
	
	
	
	done
	
	
	
	### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
	
	
	#Result folder
	mkdir Filtered_predictions	
	
	for file in Genes_predictions/*.exonerate ; do
	
		#extract some infos from file name
		file_name=`echo "$file" | sed 's/.*\///g'`
		file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
		fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
		initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
	
		#Test if the predicted gene is a Opsins gene or not
		awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
		transeq predicted_cds.fa predicted_cds.prot
		blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
		
	
		#Lets continue only if the best match is an OR
		#if grep -q -i "olfactory\|odorant" blastp_result ; then 
		if grep -q -i "RGR_\|RPE-retinal\|opsin" blastp_result ; then
	
			#Define the scaffold  name
			scaffold=`echo "$file" | sed 's/.*\///g' | sed 's/-.*//g'`
			#Define the strand on which the predicted gene is
			strand=`grep "	similarity	" $file | cut -f7`
			#Define the first position of the query on the target sequence
			first_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			#Define the last position of the query on the target sequence
			second_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			
			
			#Lets extract CDS if the gene is on the negative strand
			if [ $strand == "-" ] ; then 
			
			#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
	
	
				#If strand is minus, then the first position is:
				target_end=$((first_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((first_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
			
				#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
	
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
				query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
				perc80_query_length=$((query_length*80/100*3))
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3
				if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
			
	
	
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
			
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
					##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<5000 bp)
					if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon_rev.txt ] && rm Current_exon_rev.txt
			
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
						revseq Current_exon.fa Current_exon_rev.fa
			
						
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
			
					done
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
			
				fi
	
	
	
			#Lets make the same steps with slight modifications for the + strand
			elif [ $strand == "+" ] ; then 
			
				#If strand is minus, then the first position is:
				target_end=$((second_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((second_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
				samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
			
				#Extract the CDS sequence predicted by exonerate and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				query_name_exo=`grep -m1 "Query:" $file | sed 's/.*Query: //g'`
				query_length=`grep "$query_name_exo" Database/all_opsins.prot.fai | cut -f2`
				perc80_query_length=$((query_length*80/100*3))
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $perc80_query_length -find 3 -reverse FALSE
			
			
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
				##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" Database/all_opsins.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<5000 bp)
					if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon.txt ] && rm Current_exon.txt
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon.fa >> Current_exon.txt
			
					done
			
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
				fi
			fi
	
		else echo "$initial_header" >> No_opsins_genes_coordinates.txt
	
		fi
	
	done
fi



#Now that we have filtered all our results, we can concatenate the results

for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa

#Extract coordinates of found genes
grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
sed 's/-/	/g' No_opsins_genes_coordinates.txt >> Coordinates_already_examined.tsv
if [ `wc -l < Coordinates_already_examined.tsv` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi

current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`

#re-process blast result to find potential Opsins regions exclusing already found genes


cp Filtered_predictions/* COFFRE_PREDICTIONS/
cp Genes_predictions/* COFFRE_PREDICTIONS/


rm -r Filtered_predictions/
rm -r Genes_predictions/
rm Parsed_exonerate_gene_regions.tsv



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Filter genes and pseudogenes with a phylogenetic tree  ######################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#Remove sequences with same headers before 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_multiple_exon.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_multiple_exon_uniq.fa ; mv Pseudogenes_multiple_exon_uniq.fa Pseudogenes_multiple_exon.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshift_less_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Frameshift_less_Pseudogenes_uniq.fa ; mv Frameshift_less_Pseudogenes_uniq.fa Frameshift_less_Pseudogenes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Potential_multiple_exon_CDS.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Potential_multiple_exon_CDS_uniq.fa ; mv Potential_multiple_exon_CDS_uniq.fa Potential_multiple_exon_CDS.fa


# We now have all predicted genes in the files Pseudogenes_multiple_exon.fa and Potential_multiple_exon_CDS.fa

#First lets verify our genes with a blastp or blastxs against the uniprot database


grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' > Unverified_pseudo_id.txt

fasta_formatter -i Pseudogenes_multiple_exon.fa  > Pseudogenes_multiple_exon_reformat.fa 

transeq Potential_multiple_exon_CDS.fa Potential_multiple_exon_CDS.prot ; sed -i 's/_1//g' Potential_multiple_exon_CDS.prot
blastp -query Potential_multiple_exon_CDS.prot -db uniprot_sprot.fasta -outfmt "6 qseqid sseqid stitle" -num_threads $number_of_thread -max_target_seqs 1 -out Complete_vs_Uniprot.blastp
blastx -query Pseudogenes_multiple_exon_reformat.fa -db uniprot_sprot.fasta -outfmt "6 qseqid sseqid stitle" -num_threads $number_of_thread -max_target_seqs 1 -out Incomplete_vs_Uniprot.blastp

grep -i "RGR_\|RPE-retinal\|opsin" Complete_vs_Uniprot.blastp | cut -f1 | sort | uniq > good_complete_seq
grep -i "RGR_\|RPE-retinal\|opsin" Incomplete_vs_Uniprot.blastp | cut -f1 | sort | uniq > good_incomplete_seq

xargs samtools faidx Potential_multiple_exon_CDS.fa < good_complete_seq > Functionnal_Opsins_genes.fa
xargs samtools faidx Pseudogenes_multiple_exon_reformat.fa < good_incomplete_seq > Pseudogenes_Opsins_genes.fa


#Now lets do a rough classification of our genes using a blast

transeq Functionnal_Opsins_genes.fa Functionnal_Opsins_genes.prot ; sed -i 's/_1$//g' Functionnal_Opsins_genes.prot 
blastp -query Functionnal_Opsins_genes.prot -db Database/all_opsins.prot -outfmt "6 qseqid sseqid" -num_threads $number_of_thread -max_target_seqs 1 -out Complete_vs_Opsins.blastp
blastx -query Pseudogenes_Opsins_genes.fa -db Database/all_opsins.prot -outfmt "6 qseqid sseqid" -num_threads $number_of_thread -max_target_seqs 1 -out Incomplete_vs_Opsins.blastp

grep ">" Functionnal_Opsins_genes.fa | sed 's/>//g' | sort | uniq > uniq_id
for i in `cat uniq_id` ; do
	class=`grep -m1 "$i" Complete_vs_Opsins.blastp | cut -f2 | sed 's/-.*//g'`
	sed -i "s/>$i/>$class-$i/g" Functionnal_Opsins_genes.fa 
done



grep ">" Pseudogenes_Opsins_genes.fa | sed 's/>//g' | sort | uniq > uniq_id
for i in `cat uniq_id` ; do
	class=`grep -m1 "$i" Incomplete_vs_Opsins.blastp | cut -f2 | sed 's/-.*//g'`
	sed -i "s/>$i/>$class-$i/g" Pseudogenes_Opsins_genes.fa 
done





## Now lets filter these files : remove ambigous sequences 


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_Opsins_genes.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Functionnal_Opsins_genes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_Opsins_genes.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Pseudogenes_Opsins_genes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_Opsins_genes.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Functionnal_Opsins_genes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_Opsins_genes.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Pseudogenes_Opsins_genes.fa
cat unclear_Functionnal_Opsins_genes.fa unclear_Pseudogenes_Opsins_genes.fa > FINAL_Ambigous_Opsins.fa




## Now lets keep only one sequece of some genes are overlapping (this can happen to the extensions before predictions, but should be rare) ! 


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' FINAL_Ambigous_Opsins.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Ambiguous.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Pseudogenes_Opsins_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Functionnal_Opsins_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Functionnals.fa


nb_seq=`grep -c ">" Ambiguous.fa` ; if [ "$nb_seq" -gt "0" ] ; then grep ">" Ambiguous.fa | sed 's/>//g' | sed 's/-/,/g' | cut -f2,3,4 -d ',' > Coordinates_ambigous_final.tsv ; fi
nb_seq=`grep -c ">" Pseudogenes.fa` ; if [ "$nb_seq" -gt "0" ] ; then grep ">" Pseudogenes.fa | sed 's/>//g' | sed 's/-/,/g' | cut -f2,3,4 -d ',' > Coordinates_pseudogenes_final.tsv ; fi
nb_seq=`grep -c ">" Functionnals.fa` ; if [ "$nb_seq" -gt "0" ] ; then grep ">" Functionnals.fa | sed 's/>//g' | sed 's/-/,/g' | cut -f2,3,4 -d ',' > Coordinates_genes_final.tsv ; fi


Rscript $scripts_location/Remove_redundancy_final_dataset.R



IFS=$'\n'


for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Functionnals.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Ambiguous.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Pseudogenes.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done


if test -f "functionnal_to_keep.txt"  ; then xargs samtools faidx Functionnals.fa < functionnal_to_keep.txt > FINAL_Functional_Opsins.fa ; else echo "" > FINAL_Functional_Opsins.fa ; fi
if test -f "ambigous_to_keep.txt" ; then xargs samtools faidx Ambiguous.fa < ambigous_to_keep.txt > FINAL_Ambiguous_Opsins.fa ; else echo "" > FINAL_Ambiguous_Opsins.fa ; fi 
if test -f "pseudogenes_to_keep.txt" ; then xargs samtools faidx Pseudogenes.fa < pseudogenes_to_keep.txt > FINAL_Pseudogene_Opsins.fa ; else echo "" > FINAL_Pseudogene_Opsins.fa ; fi 


#FINAL RESULTS FILES:
#FINAL_Functional_Opsins.fa
#FINAL_Pseudogene_Opsins.fa
#FINAL_Ambiguous_Opsins.fa


nb_functionnal=`grep -c ">" FINAL_Functional_Opsins.fa`
nb_pseudo_edge=`grep -c ">" FINAL_Pseudogene_Opsins.fa`
nb_ambigous=`grep -c ">" FINAL_Ambiguous_Opsins.fa`

echo "Search of Opsins is finished. There are $nb_functionnal potentially functionnal genes, $nb_pseudo_edge pseudogenes or fragments and $nb_ambigous ambigous sequences"

echo "$nb_functionnal	$nb_pseudo_edge	$nb_ambigous" > Results_NbF_NbP_NbA_summary.txt


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"
