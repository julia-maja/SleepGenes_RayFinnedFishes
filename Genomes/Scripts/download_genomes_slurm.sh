#!/bin/bash
# Job name:
SBATCH --job-name=download_genomes
#
# Account:
SBATCH --account=def-mshafer
#
# Partition:
#SBATCH --partition=partition_name
#
# Request one node:
#SBATCH --nodes=1
#
# Specify number of tasks for use case (example):
#SBATCH --ntasks-per-node=20
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
SBATCH --time=01:00:30

## Command(s) to run (example):

sh ../index_genomes.sh
