#!/bin/bash

#SBATCH --job-name=opsin_miner 
#SBATCH --account=def-mshafer
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPUs per task
#SBATCH --mem-per-cpu=4G
#SBATCH --qos=2days
#SBATCH --time=06:00:00               # Time limit hrs:min:sec
#SBATCH --output=job_%opsin_miner.out           # Standard output file
#SBATCH --error=job_%jopsin_miner.err            # Standard error file

# Load necessary modules

# Run your command herei
	 ./Opsin_miner.sh GCF_000002035.6_GRCz11_genomic.fna one_opsin.prot uniprot_10.fasta Database/Scripts_opsins/ 50000 4
