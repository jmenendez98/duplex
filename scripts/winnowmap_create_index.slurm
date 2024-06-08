#!/bin/bash
#SBATCH --job-name=winnnowmap.index
#SBATCH --partition=long
#SBATCH --mail-user=jmmenend@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=250gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=winnowmap.index.%j.log
#SBATCH --time=168:00:00

ref_fasta=$1
ref_name="$(basename -- $ref_fasta)"
ref_name="${ref_name%.*}"

source /private/home/jmmenend/software/anaconda3/etc/profile.d/conda.sh

conda activate /private/groups/migalab/jmmenend/.conda_envs/winnowmap
meryl count k=15 output merylDB_${ref_name}_k15 ${ref_fasta}
meryl print greater-than distinct=0.9998 merylDB_${ref_name}_k15 > ${ref_name}_repetitive_k15.txt
conda deactivate

echo "Index Created: ${ref_name}_repetitive_k15.txt"