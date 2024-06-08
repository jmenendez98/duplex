#!/bin/bash
#SBATCH --job-name=pod5.convert
#SBATCH --partition=medium
#SBATCH --mail-user=jmmenend@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=pod5.convert.%j.log
#SBATCH --time=08:00:00

fast5_dir=$1

# move fast5s to a folder
mkdir -p ${fast5_dir}_fast5s
find ${fast5_dir} -name *.fast5 -exec mv {} ${fast5_dir}_fast5s \;

# pod5 convert fast5
mkdir -p ${fast5_dir}_pod5
pod5 convert fast5 ${fast5_dir}_fast5s/*.fast5 -r -o ${fast5_dir}_pod5/${fast5_dir}.pod5
