#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Download all annotated ray-finned fishes genomes on NCBI ##########################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

set -x  # Start debugging

##Date = 9 Mars 2024
##NCBI tax id for ray-finned fishes : 7898
#git clone https://github.com/pirovc/genome_updater.git
#./genome_updater/genome_updater.sh -d "refseq,genbank" -T "7898" -f "assembly_report.txt" -t 1 -m -A 1 -o Actinopterygii_assembly_reports -L curl ### HERE STOP THE COMMAND BEFORE GENOME UPDATER DOWNLOAD ALL THE GENOMES ITSELF. 
#cp Actinopterygii_assembly_reports/assembly_summary.txt ./assembly_summary.txt
#
#
##Remove partial assemblies 
#awk -v FS="\t" '$14=="Full"' assembly_summary.txt > full_assembly_summary.txt
#
#
### Now lets download all the genomes, and the annotation if it exist
#
#
#mkdir Genomic_data
cd Genomic_data
#cp ../full_assembly_summary.txt ./
#
#
module purge 
module load StdEnv/2020
module load samtools/1.15.1

IFS=$'\n'
for line in `cat full_assembly_summary_1.txt` ; do

	species=`echo "$line" | cut -f8 | sed "s/'//g" | sed 's/\.//g' | sed 's/ /_/g'`
	
	echo "$line" > line.txt
	
	rm -r $species ; mkdir $species
	
	if grep -q "GCF_" line.txt ; then 
	
		GCA_name=`echo "$line" | cut -f1`
		GCF_name=`echo "$line" | cut -f18`
		assembly_link=`echo "$line" | cut -f20 | sed 's/$/\//g' | sed "s/$GCA_name/$GCF_name/g" | sed 's/\/GCA\//\/GCF\//g'`
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		GFF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gff.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		report_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_assembly_report.txt/g' | sed "s/$GCA_name/$GCF_name/g"`
		
		report_full_link=`echo "$assembly_link$report_name"`
		full_fasta_link=`echo "$assembly_link$assembly_name"`
		full_gff_link=`echo "$assembly_link$GFF_name"`
		
		
		cd $species
		
		wget $full_gff_link ; wget $full_fasta_link ; wget $report_full_link
		gzip -d $assembly_name ; gzip -d $GFF_name 



		if grep -q "alt-scaffold" $report_name ; then

			grep -v "alt-scaffold" $report_name | cut -f7 | grep -v "^na$" | grep -v "#" | grep -v "RefSeq-Accn" | sed '/^[[:space:]]*$/d' > good_scaffold_list.txt
			grep  "alt-scaffold" $report_name | grep -v "^na$" > infos_removed_genome_parts.txt

			trimmed_genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.trimmed.fna/g'`
			genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.fna/g'`
#			xargs samtools faidx $genome_name < good_scaffold_list.txt > $trimmed_genome_name
		#	gzip $genome_name


			echo "$genome_name" >> ../List_of_modified_genomes.txt

		fi


		
		cd ../
	
	
	else
	
		assembly_link=`echo "$line" | cut -f20 | sed 's/$/\//g'`
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna.gz/g'`
		GFF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gff.gz/g'`
		report_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_assembly_report.txt/g'`

		
		full_fasta_link=`echo "$assembly_link$assembly_name"`
		full_gff_link=`echo "$assembly_link$GFF_name"`
		report_full_link=`echo "$assembly_link$report_name"`

		
		cd $species
		

		wget $full_gff_link ; wget $full_fasta_link ; wget $report_full_link
		gzip -d $assembly_name ; gzip -d $GFF_name 
		

		if grep -q "alt-scaffold" $report_name ; then

	
			grep -v "alt-scaffold" $report_name | cut -f5 | grep -v "^na$" | grep -v "#" | grep -v "GenBank-Accn" | sed '/^[[:space:]]*$/d' > good_scaffold_list.txt
			grep  "alt-scaffold" $report_name | grep -v "^na$" > infos_removed_genome_parts.txt

			trimmed_genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.trimmed.fna/g'`
			genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.fna/g'`
#			xargs samtools faidx $genome_name < good_scaffold_list.txt > $trimmed_genome_name
#			gzip $genome_name

			echo "$genome_name" >> ../List_of_modified_genomes.txt



		fi


		cd ../
	
	fi


	reduced_assembly_name=`echo "$assembly_name" | sed 's/_genomic.fna.gz//g'`
	if ls -l $species/ | grep -q ".gff" ; then 
		echo "$species,$reduced_assembly_name" >> Annotated_genome_list.txt
	else 
		echo "$species,$reduced_assembly_name" >> Non_annotated_genome_list.txt
	fi 


done



