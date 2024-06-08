#!/bin/bash
#SBATCH --job-name=winnnowmap.ONT
#SBATCH --partition=long
#SBATCH --mail-user=jmmenend@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=250gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=winnowmap.ONT.%j.log
#SBATCH --time=168:00:00

set -eux -o pipefail

source /private/home/jmmenend/software/anaconda3/etc/profile.d/conda.sh

methyl_bam=$1
output_prefix="$(basename ${methyl_bam%.bam})"

ref_fasta=$2
ref_name="$(basename -- $ref_fasta)"
ref_name="${ref_name%.*}"

temp="/private/groups/migalab/jmmenend/tmp/"

methyl_fastq="${output_prefix}.fastq"
samtools fastq -@ 64 -T MM,ML,Mm,Ml ${methyl_bam} > ${methyl_fastq}

conda activate /private/groups/migalab/jmmenend/.conda_envs/winnowmap
#meryl count k=15 output merylDB_${ref_name}_k15 ${ref_fasta}
#meryl print greater-than distinct=0.9998 merylDB_${ref_name}_k15 > ${ref_name}_repetitive_k15.txt
winnowmap \
	--cs \
	--MD \
	-I 16G \
	-W ${ref_name}_repetitive_k15.txt \
	-y \
	-k 15 \
	-ax map-ont \
	${ref_fasta} ${methyl_fastq} \
	> ${output_prefix}_winnowmap.sam
conda deactivate

samtools view -bh -@ 64 ${output_prefix}_winnowmap.sam | \
	samtools sort -@ 64 - > ${output_prefix}_winnowmap_sort.bam
samtools index -@ 64 ${output_prefix}_winnowmap_sort.bam

rm ${output_prefix}_winnowmap.sam
rm ${methyl_fastq}

echo "Done."

