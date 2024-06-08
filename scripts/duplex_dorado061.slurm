#!/bin/bash
#SBATCH --job-name=dorado.duplex
#SBATCH --partition=gpu
#SBATCH --mail-user=jmmenend@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=1000gb
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=64
#SBATCH --output=dorado.duplex.%j.log
#SBATCH --time=168:00:00

pod5_folder=$1
output=$2

/private/home/jmmenend/bin/dorado-0.6.1 duplex sup,5mCG_5hmCG $pod5_folder > $output

