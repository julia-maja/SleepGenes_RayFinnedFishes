##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Clean Reads  ###################################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

mkdir FastQC_Results ; mkdir FastP_results

cut -f2 RNA_info_table.tsv > list_sp #list of species with RNA_seq data.

for species in `cat list_sp` ; do

	cd $species

	cp ../Fastp_RNA_reads.sh ./

	ls -l | grep "fastq" | sed 's/.* //g' | sed 's/.R1.*//g' | sed 's/.R2.*//g' | tail -n+2 | sort | uniq > list_tissues_reads

	for tissue_reads in `cat list_tissues_reads` ; do

		sbatch --qos=6hours -c 4 --mem=10G Fastp_RNA_reads.sh $tissue_reads

	done

	cd ../

done





========================= Fastp_RNA_reads.sh  ===========================================================================

#!/bin/bash


#SBATCH --job-name=Fastp_RNA_reads   # Job name


eval "$(conda shell.bash hook)"
conda activate olfactory


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load fastp #(https://onlinelibrary.wiley.com/doi/10.1002/imt2.107)
module load FastQC


tissue=$1 


R1_reads=`echo "$tissue.R1.fastq.gz"`
R2_reads=`echo "$tissue.R2.fastq.gz"`

R1_reads_cleaned=`echo "$R1_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`
R2_reads_cleaned=`echo "$R2_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`

report_file=`echo "$tissue.report"`


fastp -i $R1_reads -I $R2_reads -o $R1_reads_cleaned -O $R2_reads_cleaned -h $report_file.html -j $report_file.json

======================================================================================================================================================

###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
########################################## Create genome index files for STAR  ######################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#Extract genomes for which I have RNA-seq data

cp ../ALL_Genomes/full_genome_assembly_summary.txt ./
cut -f2 RNA_info_table.tsv > list_sp 

rm RNAseq_genome_assembly_summary.txt
for species in `cat list_sp` ; do 
	species_name=`echo "$species" | sed 's/_/ /g'`

	grep "	$species_name	" full_genome_assembly_summary.txt >> RNAseq_genome_assembly_summary.txt
done 


#Download FASTA, GTF and GFF3 files 

mkdir GTF_files


IFS=$'\n'
for line in `cat RNAseq_genome_assembly_summary.txt` ; do

	species=`echo "$line" | cut -f8 | sed "s/'//g" | sed 's/\.//g' | sed 's/ /_/g'`
	
	echo "$line" > line.txt
	
	rm -r $species.Genome ; mkdir $species.Genome
	
	if grep -q "GCF_" line.txt ; then 
	
		GCA_name=`echo "$line" | cut -f1`
		GCF_name=`echo "$line" | cut -f18`
		assembly_link=`echo "$line" | cut -f20 | sed 's/$/\//g' | sed "s/$GCA_name/$GCF_name/g" | sed 's/\/GCA\//\/GCF\//g'`
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		GFF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gff.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		GTF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gtf.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		report_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_assembly_report.txt/g' | sed "s/$GCA_name/$GCF_name/g"`
		
		report_full_link=`echo "$assembly_link$report_name"`
		full_fasta_link=`echo "$assembly_link$assembly_name"`
		full_gff_link=`echo "$assembly_link$GFF_name"`
		full_gtf_link=`echo "$assembly_link$GTF_name"`
		
		
		cd $species.Genome
		
		wget $full_gff_link ; wget $full_fasta_link ; wget $report_full_link ; wget $full_gtf_link 
		
		
		gzip -d $assembly_name ; gzip -d $GFF_name ; gzip -d $GTF_name 
		

		if grep -q "alt-scaffold" $report_name ; then

			grep -v "alt-scaffold" $report_name | cut -f7 | grep -v "^na$" | grep -v "#" | grep -v "RefSeq-Accn" | sed '/^[[:space:]]*$/d' > good_scaffold_list.txt
			grep  "alt-scaffold" $report_name | grep -v "^na$" > infos_removed_genome_parts.txt

			trimmed_genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.trimmed.fna/g'`
			genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.fna/g'`
			xargs samtools faidx $genome_name < good_scaffold_list.txt > $trimmed_genome_name
			gzip $genome_name


			echo "$genome_name" >> ../List_of_modified_genomes.txt

		fi


		
		cd ../
	
	
	else
	
		assembly_link=`echo "$line" | cut -f20 | sed 's/$/\//g'`
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna.gz/g'`
		GFF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gff.gz/g'`
		GTF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gtf.gz/g'`
		report_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_assembly_report.txt/g'`

		
		full_fasta_link=`echo "$assembly_link$assembly_name"`
		full_gff_link=`echo "$assembly_link$GFF_name"`
		full_gtf_link=`echo "$assembly_link$GTF_name"`
		report_full_link=`echo "$assembly_link$report_name"`

		
		cd $species.Genome
		
		wget $full_gff_link ; wget $full_fasta_link ; wget $report_full_link ; wget $full_gtf_link
		
		
		gzip -d $assembly_name ; gzip -d $GFF_name ; gzip -d $GTF_name 


		if grep -q "alt-scaffold" $report_name ; then

	
			grep -v "alt-scaffold" $report_name | cut -f5 | grep -v "^na$" | grep -v "#" | grep -v "GenBank-Accn" | sed '/^[[:space:]]*$/d' > good_scaffold_list.txt
			grep  "alt-scaffold" $report_name | grep -v "^na$" > infos_removed_genome_parts.txt

			trimmed_genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.trimmed.fna/g'`
			genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.fna/g'`
			xargs samtools faidx $genome_name < good_scaffold_list.txt > $trimmed_genome_name
			gzip $genome_name

			echo "$genome_name" >> ../List_of_modified_genomes.txt



		fi



		
		cd ../
	
	fi


done



#Make a table with genomic info per species

IFS=$'\n'
rm RNAspecies_Info_files.csv 
for line in `cat RNAseq_genome_assembly_summary.txt` ; do

	species=`echo "$line" | cut -f8 | sed "s/'//g" | sed 's/\.//g' | sed 's/ /_/g'`


	ls -l $species.Genome | sed 's/.* //g' > current_ls.txt 

	if grep -q ".gtf$" current_ls.txt  ; then 

		GTF_file=`grep ".gtf$" current_ls.txt` 
		GTF_full_path=`echo "/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/RNA_SEQ_Analysis/$species.Genome/$GTF_file"`

	else

		GTF_full_path=""

	fi 

	if grep -q ".gff$" current_ls.txt  ; then 

		GFF_file=`grep ".gff$" current_ls.txt` 
		GFF_full_path=`echo "/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/RNA_SEQ_Analysis/$species.Genome/$GFF_file"`

	else 

		GFF_full_path=""

	fi 

	Genome_file=`grep ".fna$" current_ls.txt` 
	Genome_full_path=`echo "/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/RNA_SEQ_Analysis/$species.Genome/$Genome_file"`

	Info_file=`grep "report.txt$" current_ls.txt` 
	Info_full_path=`echo "/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/RNA_SEQ_Analysis/$species.Genome/$Info_file"`


	
	echo "$species,$Genome_full_path,$Info_full_path,$GTF_full_path,$GFF_full_path" >> RNAspecies_Info_files.csv 


done



# Now use those info to launch STAR. For species with a gtf file, we can use it on STAR to improve alignment .. 

IFS=$'\n'
for line in `cat RNAspecies_Info_files.csv` ; do 

	echo "$line" > current_line.txt 

	if grep -q "gtf," current_line.txt ; then 

		species_name=`cut -f1 -d "," current_line.txt`
		genome_dir=`echo "$species_name.Genome"`
		genome_file=`cut -f2 -d "," current_line.txt`
		GTF_file=`cut -f4 -d "," current_line.txt`


		sbatch -e error.$species_name.out -o slurm.$species_name.out --qos=1day -c 8 --mem=50G index_genome_star_annot.sh $genome_dir $genome_file $GTF_file 13

	else 


		species_name=`cut -f1 -d "," current_line.txt`
		genome_dir=`echo "$species_name.Genome"`
		genome_file=`cut -f2 -d "," current_line.txt`

		sbatch -e error.$species_name.out -o slurm.$species_name.out --qos=1day -c 8 --mem=50G index_genome_star.sh $genome_dir $genome_file 13
	fi


	rm current_line.txt 

done




###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
########################################## Create GFF3 files for each genome ########################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

#First import GFF3 from genome annotations 

IFS=$'\n'
for line in `cat RNAspecies_Info_files.csv` ; do 

	echo "$line" > current_line.txt 
	species=`cut -f1 -d "," current_line.txt `

	if grep -q "gtf," current_line.txt ; then 

		GFF_parsed_name=`ls -l /scicore/home/salzburg/fogg0000/Vertebrate_genomes/$species/ | grep "gff.LI" | sed 's/.* //g'`
		cp /scicore/home/salzburg/fogg0000/Vertebrate_genomes/$species/$GFF_parsed_name $species.Genome/$species.parsed.gff3

		echo "/scicore/home/salzburg/fogg0000/Vertebrate_genomes/$species/$GFF_parsed_name"

	fi 

done


#Now create opsins GFF3 files


IFS=$'\n'

for species in `cat list_sp` ; do

	rm -r $species.GFF3.opsins/ ; mkdir $species.GFF3.opsins/
	grep ">$species---" /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/TOTAL_opsins.fa | sed 's/>//g' | cut -f5,6,7 -d "-" | sed 's/-/,/g' > $species.GFF3.opsins/$species.opsins.csv

	genome_fasta=`grep "$species" RNAspecies_Info_files.csv | cut -f2 -d ","`
	annotation_rerort_file=`grep "$species" RNAspecies_Info_files.csv  | cut -f3 -d ","`
	
	for line in `cat $species.GFF3.opsins/$species.opsins.csv` ; do
		scaffold=`echo "$line" | cut -f1 -d ","`
		start=`echo "$line" | cut -f2 -d ","`
		stop=`echo "$line" | cut -f3 -d ","`
	
	
		if grep -q "$species.*GCF_" RNAspecies_Info_files.csv ; then
			new_scaff=`grep "$scaffold" $annotation_rerort_file | cut -f7`
		else 
			new_scaff=$scaffold
		fi
	
		start_ext=$((start - 1000))
		stop_ext=$((stop + 1000))
		if [ "$start_ext" -lt 1 ] ; then start_ext=1; fi
	
		samtools faidx $genome_fasta $new_scaff:$start_ext-$stop_ext > $species.GFF3.opsins/current_scaff.fa
	
		gene_name=`grep "$scaffold.*$start.*$stop" /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/TOTAL_opsins.fa | sed 's/>//g'`
		samtools faidx /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/TOTAL_opsins.fa $gene_name > $species.GFF3.opsins/current_gene.fa
	
	
		/scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/Custom_GFF3/exonerate-gff3/src/program/exonerate --model est2genome --bestn 1 --showtargetgff TRUE --gff3 TRUE $species.GFF3.opsins/current_gene.fa $species.GFF3.opsins/current_scaff.fa > $species.GFF3.opsins/$gene_name.exonerate
	
	
	done 
	
	
	
	for file in $species.GFF3.opsins/*.exonerate ; do 
		
		file_name=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
	
		grep "exonerate:est2genome" $file > $species.GFF3.opsins/$file_name.gff3
		sed -i 's/match_part/exon/g' $species.GFF3.opsins/$file_name.gff3
		sed -i "s/match00001/$file_name/g" $species.GFF3.opsins/$file_name.gff3
		sed -i "s/match00002/$file_name/g" $species.GFF3.opsins/$file_name.gff3
		sed -i 's/exonerate:est2genome/exonerate/g' $species.GFF3.opsins/$file_name.gff3
		sed -i 's/match/gene/g' $species.GFF3.opsins/$file_name.gff3
	
		start_scaff=`cut -f1 $species.GFF3.opsins/$file_name.gff3 | head -1 | sed 's/:/-/g' | cut -f2 -d "-"` 
		cut -f4 $species.GFF3.opsins/$file_name.gff3 > $species.GFF3.opsins/start_fatures
		cut -f5 $species.GFF3.opsins/$file_name.gff3 > $species.GFF3.opsins/stop_fatures
	
	
		rm $species.GFF3.opsins/start_fatures_true
		for start_pos in `cat $species.GFF3.opsins/start_fatures` ; do 
			start_pos_new=$((start_pos + start_scaff))
			echo "$start_pos_new" >> $species.GFF3.opsins/start_fatures_true
		done
		
		rm $species.GFF3.opsins/stop_fatures_true
		for stop_pos in `cat $species.GFF3.opsins/stop_fatures` ; do 
			stop_pos_new=$((stop_pos + start_scaff))
			echo "$stop_pos_new" >> $species.GFF3.opsins/stop_fatures_true
		done
		
		cut -f1 $species.GFF3.opsins/$file_name.gff3 | sed 's/:.*//g' > $species.GFF3.opsins/col1
		cut -f2,3 $species.GFF3.opsins/$file_name.gff3 | sed 's/:.*//g' > $species.GFF3.opsins/col2-3
		cut -f6- $species.GFF3.opsins/$file_name.gff3 > $species.GFF3.opsins/col6
		paste -d "\t" $species.GFF3.opsins/col1 $species.GFF3.opsins/col2-3 > $species.GFF3.opsins/col1-3
		paste -d "\t" $species.GFF3.opsins/col1-3 $species.GFF3.opsins/start_fatures_true > $species.GFF3.opsins/col1-4
		paste -d "\t" $species.GFF3.opsins/col1-4 $species.GFF3.opsins/stop_fatures_true > $species.GFF3.opsins/col1-5
		paste -d "\t" $species.GFF3.opsins/col1-5 $species.GFF3.opsins/col6 > $species.GFF3.opsins/$file_name.gff3.parsed
	
	
		gene_line=`head -1 $species.GFF3.opsins/$file_name.gff3.parsed`
		echo "$gene_line" > $species.GFF3.opsins/gene_line.txt
		sed -i 's/;.*//g' $species.GFF3.opsins/gene_line.txt
		gene_ID=`cut -f9 $species.GFF3.opsins/gene_line.txt | sed 's/ID=//g'`
		sed 's/	gene	/	mRNA	/g' $species.GFF3.opsins/gene_line.txt > $species.GFF3.opsins/RNA_line.txt
		sed -i 's/ID=/ID=RNA-/g' $species.GFF3.opsins/RNA_line.txt
		sed -i "s/$/;Parent=$gene_ID;Name=TAAR-$gene_ID;/g" $species.GFF3.opsins/RNA_line.txt
		sed -i 's/$/;/g' $species.GFF3.opsins/gene_line.txt
	
		grep "	exon	" $species.GFF3.opsins/$file_name.gff3.parsed > $species.GFF3.opsins/exon_lines.txt
	
	
		i=1
		rm $species.GFF3.opsins/exon_lines_parsed.txt
		for exon_line in `cat $species.GFF3.opsins/exon_lines.txt` ; do 
			echo "$exon_line" | sed "s/Parent=.*/ID=exon-$gene_ID-$i;Parent=RNA-$gene_ID;/g" >> $species.GFF3.opsins/exon_lines_parsed.txt
			i=$((i+1))
		done
		
		cp $species.GFF3.opsins/exon_lines_parsed.txt $species.GFF3.opsins/CDS_lines_parsed.txt
		sed -i 's/	exon	/	CDS	/g' $species.GFF3.opsins/CDS_lines_parsed.txt
		
		
		
		sed -i "s/.	ID=exon-$gene_ID.*;Parent=/.	ID=cds-$gene_ID;Parent=/g" $species.GFF3.opsins/CDS_lines_parsed.txt
		
		cat $species.GFF3.opsins/gene_line.txt $species.GFF3.opsins/RNA_line.txt $species.GFF3.opsins/exon_lines_parsed.txt $species.GFF3.opsins/CDS_lines_parsed.txt > $species.GFF3.opsins/$file_name.gff3.parsed
	
	
	done

	#accession=`grep "$species" RNAspecies_Info_files.csv | cut -f2 -d "," | sed 's/.*Genome.//g' | sed 's/GCA_/GCA-/g' | sed 's/GCF_/GCF-/g' | cut -f1 -d "_" | sed 's/-/_/g'`
	#name_build=`grep "$species" RNAspecies_Info_files.csv | cut -f2 -d "," | sed 's/.*Genome.//g' | sed 's/GCA_/GCA-/g' | sed 's/GCF_/GCF-/g' | cut -f2 -d "_"`

	echo "##gff-version 3" > $species.GFF3.opsins/opsins.gff3
	echo '#!gff-spec-version 1.21' >> $species.GFF3.opsins/opsins.gff3
	echo '#!processor NCBI annotwriter' >> $species.GFF3.opsins/opsins.gff3
	echo '#!genome-build current_build' >> $species.GFF3.opsins/opsins.gff3
	echo '#!genome-build-accession NCBI_Assembly:accession' >> $species.GFF3.opsins/opsins.gff3
	echo '#!annotation-source NCBI name_build' >> $species.GFF3.opsins/opsins.gff3
	echo '##sequence-region scaffold 1 40673430' >> $species.GFF3.opsins/opsins.gff3
	echo '##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=8128' >> $species.GFF3.opsins/opsins.gff3
	cat $species.GFF3.opsins/*.gff3.parsed >> $species.GFF3.opsins/opsins.gff3
	sed -i 's/=exon-/=opsins-exon-/g' $species.GFF3.opsins/opsins.gff3 


done



IFS=$'\n'

for species in `cat list_sp` ; do

	x=`grep -c "	gene	" $species.GFF3.opsins/opsins.gff3`
	y=`grep -c ">$species---" /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Phylogenies_Opsins_v2/TOTAL_opsins.fa`

	echo "$species,$x,$y"

done


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
########################################## Align reads to genomes   ##############################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

#now for all species


for species in `cat list_sp` ; do

	cd $species

	cp ../Align_STAR.sh ./

	rm slurm*

	ls -l | grep "fastq" | sed 's/.* //g' | sed 's/.R1.*//g' | sed 's/.R2.*//g' | tail -n+2 | sort | uniq > list_tissues_reads

	for tissue in `cat list_tissues_reads` ; do

		R1_reads=`echo "$tissue.R1.fastq.gz"`
		R2_reads=`echo "$tissue.R2.fastq.gz"`
		R1_reads_cleaned=`echo "$R1_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`
		R2_reads_cleaned=`echo "$R2_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`


		sbatch --qos=6hours -c 40 --mem=150G --job-name=$tissue Align_STAR.sh $species $tissue $R1_reads_cleaned $R2_reads_cleaned

	done 


	cd ../

done




for species in `cat list_sp` ; do ls -lh $species | grep ".bam$" ; done 



## Verify that I have a bam file for each tissues, and that the good command was launched by STAR .. 



for species in `cat list_sp` ; do

	cd $species

	for tissue in `cat list_tissues_reads` ; do

		log_file=`echo "$tissue,Log.out" | sed 's/,//g'`

		grep -A1 "##### Command Line:" $log_file | tail -1 | sed 's/.*readFilesIn//g'  | sed 's/^ //g' | sed 's/ ..outSAMtype.*//g' | tr ' ' '\n' > reads_used
		nb_good=`grep -c "$tissue" reads_used`


		if [ $nb_good -lt 2 ] ; then echo "$tissue" ; fi

		if ls -l | grep -q "$tissue.*.bam$" ; then echo "$tissue" >> ../ok_file ; else echo "$tissue" ; fi

	done

	cd ../ 

done




## Compute stats for each run, such as the nb of reads, nb of mapped reads etc ....


rm STAR_all_statistics.csv 
for species in `cat list_sp` ; do

	cd $species 

	for tissue in `cat list_tissues_reads` ; do

		log_file=`echo "$tissue,Log.final.out" | sed 's/,//g'`
		tissue_uniq_name=`echo "$tissue" | sed 's/.*\.//g'`

		sed 's/^ *//g' $log_file | sed 's/://g' | sed 's/%/perc_/g' | sed 's/\//_/g' | sed 's/perc_$//g' | grep -v "Million of reads per hour" | grep -v "Started job on" | grep -v "Finished on" | grep -v "Started mapping on" | grep -v "CHIMERIC READS" | grep -v "MULTI-MAPPING READS" | grep -v "UNMAPPED READS" | grep -v "UNIQUE READS" | sed 's/|/,/g' | sed 's/ ,/,/g' | sed 's/,	/,/g' | sed 's/ /_/g' | sed 's/_,/,/g' | sed '/^$/d' | sed "s/^/$species,$tissue_uniq_name,/g" >> ../STAR_all_statistics.csv

	done

	cd ../

done


# Make a nice table 


IFS=$'\n'
rm  Summary_NCBI_SRA_data.csv
for line in `cat table_misc` ; do 
	order=`echo "$line" | cut -f1`
	species=`echo "$line" | cut -f2`

	eye=`echo "$line" | cut -f3 | tr ',' "\n"` ; liver=`echo "$line" | cut -f4 | tr ',' "\n"` ; intestine=`echo "$line" | cut -f5 | tr ',' "\n"` ; fin=`echo "$line" | cut -f6 | tr ',' "\n"` ; heart=`echo "$line" | cut -f7 | tr ',' "\n"` ; skin=`echo "$line" | cut -f8 | tr ',' "\n"` ; gill=`echo "$line" | cut -f9 | tr ',' "\n"` ; muscle=`echo "$line" | cut -f10 | tr ',' "\n"` ; brain=`echo "$line" | cut -f11 | tr ',' "\n"` ; kidney=`echo "$line" | cut -f12 | tr ',' "\n"` ; bones=`echo "$line" | cut -f13 | tr ',' "\n"` ; spleen=`echo "$line" | cut -f14 | tr ',' "\n"` ; testis=`echo "$line" | cut -f15 | tr ',' "\n"` ; ovary=`echo "$line" | cut -f16 | tr ',' "\n"`

	echo "$eye" | sed "s/^/$order,$species,eye,/g" >> Summary_NCBI_SRA_data.csv
	echo "$liver" | sed "s/^/$order,$species,liver,/g" >> Summary_NCBI_SRA_data.csv
	echo "$intestine" | sed "s/^/$order,$species,intestine,/g" >> Summary_NCBI_SRA_data.csv
	echo "$fin" | sed "s/^/$order,$species,fin,/g" >> Summary_NCBI_SRA_data.csv
	echo "$heart" | sed "s/^/$order,$species,heart,/g" >> Summary_NCBI_SRA_data.csv
	echo "$skin" | sed "s/^/$order,$species,skin,/g" >> Summary_NCBI_SRA_data.csv
	echo "$gill" | sed "s/^/$order,$species,gill,/g" >> Summary_NCBI_SRA_data.csv
	echo "$muscle" | sed "s/^/$order,$species,muscle,/g" >> Summary_NCBI_SRA_data.csv
	echo "$brain" | sed "s/^/$order,$species,brain,/g" >> Summary_NCBI_SRA_data.csv
	echo "$kidney" | sed "s/^/$order,$species,kidney,/g" >> Summary_NCBI_SRA_data.csv
	echo "$bones" | sed "s/^/$order,$species,bones,/g" >> Summary_NCBI_SRA_data.csv
	echo "$spleen" | sed "s/^/$order,$species,spleen,/g" >> Summary_NCBI_SRA_data.csv 
	echo "$testis" | sed "s/^/$order,$species,testis,/g" >> Summary_NCBI_SRA_data.csv 
	echo "$ovary" | sed "s/^/$order,$species,ovary,/g" >> Summary_NCBI_SRA_data.csv 


done

grep -v ",$" Summary_NCBI_SRA_data.csv > temp.csv ; mv temp.csv Summary_NCBI_SRA_data.csv



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
########################################## Create index for all bam files  ##############################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#Create an index for bam files

for species in `cat list_sp` ; do

	cd $species ; cp ../samtools_index.sh ./ ; rm slurm*

	for file in *.bam ; do sbatch --qos=6hours -c 8 --mem=8G samtools_index.sh $file ; done


	cd ../

done



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
############ Count the reads manually annotated opsin genes ##########################################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#### First on OPSINS


mkdir HTseq_Count_Results

for species in `cat list_sp` ; do

	cd $species ; cp ../HTSeq_count.sh ./

	for i in *.bam ; do 

		sbatch --qos=6hours -c 4 --mem=8G -e error.$species.out -o slurm.$species.out HTSeq_count.sh $species $i

	done



	cd ../

done



rm ALL_species_tissues_RNA.csv
cd HTseq_Count_Results

for i in *.count ; do 

	
	species=`echo "$i" | cut -f1 -d "."`
	tissue=`echo "$i" | cut -f2 -d "."`

	sed 's/opsins-exon-//g' $i | sed 's/exons-/exons,/g' | sed "s/^/$species,$tissue,/g" | sed 's/	/,/g' | head -n -5 >> ../ALL_species_tissues_RNA.csv

done




###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
############ Count the reads on automatically annotated genes ########################################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#### Now on ALL annotations (for species which have annot ....)


mkdir HTseq_Count_Results_ALL_GENES


for species in `cat list_sp` ; do grep "$species" RNAspecies_Info_files.csv >> temp.csv ; done
mv temp.csv RNAspecies_Info_files_up.csv


IFS=$'\n'
for line in `cat RNAspecies_Info_files_up.csv` ; do 

	echo "$line" > current_line.txt 
	species=`cut -f1 -d "," current_line.txt`

	if grep -q "gtf," current_line.txt ; then 

		cd $species ; cp ../HTSeq_count_allgenes.sh ./ ; rm slurm*

		for i in *.bam ; do 

			sbatch --qos=1day -c 8 --mem=50G HTSeq_count_allgenes.sh $species $i ../$species.Genome/$species.parsed.gff3

		done

		cd ../

	fi 


done






wc -l *.count | grep " 0 "




rm Annotated_species_tissues_all_genes_RNA.csv
cd HTseq_Count_Results_ALL_GENES

for i in *.count ; do 

	
	species=`echo "$i" | cut -f1 -d "."`
	tissue=`echo "$i" | cut -f2 -d "."`

	sed 's/exon-//g' $i | sed 's/-[0-9]*	/,/g' | sed 's/-[0-9]*	//g' | sed "s/^/$species,$tissue,/g"  |  sed 's/gnl.*|//g' | head -n -5 >> ../Annotated_species_tissues_all_genes_RNA.csv

done

sed -i 's/	/,/g' Annotated_species_tissues_all_genes_RNA.csv


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
############ Extarct genes corresponding to opsins + extract the OGG of each gene ###################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================



#Make a table with the species and opsin genes annotated in the automatic annotation


rm Species_AnnotatedOpsins.csv

IFS=$'\n'
for line in `cat RNAspecies_Info_files_up.csv` ; do 

	echo "$line" > current_line.txt 
	species=`cut -f1 -d "," current_line.txt`

	if grep -q "gtf," current_line.txt ; then 

		cut -f1 -d "," /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/ALL_Genomes/Automatic_annotations/$species.Sequences_Comp/annotated_locations.csv | sed "s/^/$species,/g" | sed 's/rna-//g' >> Species_AnnotatedOpsins.csv

	fi

done




#Make a table with each gene and their OGG ... 


rm Table_corresp_OGG_gene.csv
cut -f1 ~/Non_visual_opsins_Project/OGG_analysis/Proteomes/OrthoFinder/Results_Dec24/Orthogroups/Orthogroups.tsv | tail -n+2 > list_ogg

IFS=$'\n'
for line in `cat RNAspecies_Info_files_up.csv` ; do 

	echo "$line" > current_line.txt 
	species=`cut -f1 -d "," current_line.txt`

	if grep -q "gtf," current_line.txt ; then 

		echo "$species"
		awk -F"\t" -v column_val="$species" '{ if (NR==1) {val=-1; for(i=1;i<=NF;i++) { if ($i == column_val) {val=i;}}} if(val != -1) print $val} ' ~/Non_visual_opsins_Project/OGG_analysis/Proteomes/OrthoFinder/Results_Dec24/Orthogroups/Orthogroups.tsv | tail -n+2 > curr_species_col.$species

		nb_OGG=`wc -l < list_ogg`

		for index in `seq 1 $nb_OGG` ; do

			curr_OGG=`head -$index list_ogg | tail -1`
			head -$index curr_species_col.$species | tail -1 | tr ',' "\n" | sed 's/^ *//g' | sed "s/^/$curr_OGG,/g" | sed "s/$species---//g" | sed "s/rna-//g" | sed "s/_1$//g" | sed "s/^/$species,/g" >> Table_corresp_OGG_gene.csv



		done

		rm curr_genes.$species
		rm curr_species_col.$species


	fi



done




cat *.OGG.GeneName.csv > Table_corresp_OGG_gene.csv




#Make a table with the gene length per gene ... 


rm Table_all_gene_length.csv

IFS=$'\n'
for line in `cat RNAspecies_Info_files_up.csv` ; do 

	echo "$line" > current_line.txt 
	species=`cut -f1 -d "," current_line.txt`

	if grep -q "gtf," current_line.txt ; then 

		echo "$species"

		samtools faidx ../OGG_analysis/Proteomes/$species.fa

		cut -f1,2 ../OGG_analysis/Proteomes/$species.fa.fai | sed "s/$species---//g" | sed "s/rna-//g" | sed 's/_1	/,/g' | tr "\t" ',' | sed "s/^/$species,/g" >> Table_all_gene_length.csv


	fi 

done




=========================== 


#!/bin/bash


#SBATCH --job-name=Extract_OGG   # Job name

species=$1 



IFS=$'\n'

awk -F"\t" -v column_val="$species" '{ if (NR==1) {val=-1; for(i=1;i<=NF;i++) { if ($i == column_val) {val=i;}}} if(val != -1) print $val} ' OrthoFinder/Results_Dec24/Orthogroups/Orthogroups.tsv > $species.curr_species_col


if [ $species == "Xyrauchen_texanus" ] ; then cut -f179 OrthoFinder/Results_Dec24/Orthogroups/Orthogroups.tsv > $species.curr_species_col ; fi

paste -d "\t" OGG_col $species.curr_species_col | tail -n +2 > $species.curr_species_OGG


for line in `cat $species.curr_species_OGG` ; do 
	curr_OGG=`echo "$line" | cut -f1`
	echo "$line" | cut -f2 | sed 's/, /,/g' | tr , '\n' > $species.curr_species_genes


	for gene in `cat $species.curr_species_genes` ; do 

		grep "^$gene	" out_interpro/$species.fa.tsv | cut -f12 | grep -v "^-$" | sort | uniq > IPR_list.$species
		grep "^$gene	" out_interpro/$species.fa.tsv | cut -f14 | grep -v "^-$" | tr "|" '\n' | sort | uniq > Goterm_list.$species


		for IPR in `cat IPR_list.$species` ; do 
			echo "$curr_OGG,$species,$gene,$IPR" >> OGG_IPR_table.csv
		done

		for GOterm in `cat Goterm_list.$species` ; do 
			echo "$curr_OGG,$species,$gene,$GOterm" >> OGG_GO_table.csv
		done



	done
done


===========================



########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
###### ACCESORRY SCRIPT ######
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

================ samtools_index.sh ================================================================================================================

#!/bin/bash


#SBATCH --job-name=samtools_index   # Job name


eval "$(conda shell.bash hook)"
conda activate olfactory

module load SAMtools

file=$1

samtools index $file




================================================================================================================================================================






================ index_genome_star_annot.sh ================================================================================================================
#!/bin/bash


#SBATCH --job-name=STAR_index_annot   # Job name


eval "$(conda shell.bash hook)"
conda activate miniprot


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load SAMtools
module load STAR

genome_dir=$1
genome_file=$2
GTF_file=$3
indexNBase=$4

samtools faidx $genome_file

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $genome_file --sjdbGTFfile $GTF_file --genomeSAindexNbases $indexNBase


================================================================================================================================================================================



================ index_genome_star.sh ================================================================================================================
#!/bin/bash


#SBATCH --job-name=STAR_index   # Job name


eval "$(conda shell.bash hook)"
conda activate miniprot


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load SAMtools
module load STAR

genome_dir=$1
genome_file=$2
indexNBase=$3

samtools faidx $genome_file

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $genome_file --genomeSAindexNbases $indexNBase


================================================================================================================================================================================





===== Align_STAR.sh ==================================================



#!/bin/bash


eval "$(conda shell.bash hook)"
conda activate olfactory

module load STAR


species=$1
tissue=$2
R1_reads_cleaned=$3
R2_reads_cleaned=$4


STAR --runMode alignReads --runThreadN 20 --genomeDir /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/RNA_SEQ_Analysis/$species.Genome/ --readFilesIn $R1_reads_cleaned $R2_reads_cleaned --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix $tissue --outSAMunmapped Within --outSAMmultNmax 1




===============================================================




================================HTSeq_count.sh=========================================================

#!/bin/bash


#SBATCH --job-name=HTSeq_count   # Job name


eval "$(conda shell.bash hook)"
conda activate HTSEQ_env

species=$1
bam_file=$2
sample_name=`echo "$bam_file" | sed 's/Aligned.sortedByCoord.out.bam//g'`

htseq-count -f bam -r pos -m union -s reverse -t exon $bam_file /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/RNA_SEQ_Analysis/$species.GFF3.opsins/opsins.gff3 --idattr=ID > ../HTseq_Count_Results/$sample_name.opsins.count


=========================HTSeq_count_allgenes.sh===========================================================


#!/bin/bash


#SBATCH --job-name=HTSeq_count   # Job name


eval "$(conda shell.bash hook)"
conda activate HTSEQ_env

species=$1
bam_file=$2
annot_file=$3
sample_name=`echo "$bam_file" | sed 's/Aligned.sortedByCoord.out.bam//g'`

htseq-count -f bam -r pos -m union -s reverse -t exon $bam_file $annot_file --idattr=ID > ../HTseq_Count_Results_ALL_GENES/$sample_name.ALL.count



============================================================================================================





extract_ogg_ganename.sh ======== 



#!/bin/bash


#SBATCH --job-name=OGG_to_geneName   # Job name



eval "$(conda shell.bash hook)"
conda activate miniprot

species=$1

rm curr_species_col.$species
rm $species.OGG.GeneName.csv
awk -F"\t" -v column_val="$species" '{ if (NR==1) {val=-1; for(i=1;i<=NF;i++) { if ($i == column_val) {val=i;}}} if(val != -1) print $val} ' ~/Non_visual_opsins_Project/OGG_analysis/Proteomes/OrthoFinder/Results_Dec24/Orthogroups/Orthogroups.tsv | tail -n+2 > curr_species_col.$species

nb_OGG=`wc -l < list_ogg`

for index in `seq 1 $nb_OGG` ; do

	curr_OGG=`head -$index list_ogg | tail -1`
	head -$index curr_species_col.$species | tail -1 | tr ',' "\n" | sed 's/^ *//g' | sed "s/^/$curr_OGG,/g" | sed "s/$species---//g" | sed "s/rna-//g" | sed "s/_1$//g" | sed "s/^/$species,/g" >> $species.OGG.GeneName.csv



done

rm curr_species_col.$species



