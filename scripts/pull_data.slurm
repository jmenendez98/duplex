#!/bin/bash
#SBATCH --job-name=pull_duplex
#SBATCH --partition=long
#SBATCH --mail-user=jmmenend@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=pull_duplex.%j.log
#SBATCH --time=168:00:00

links=(
    "s3://human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_Duplex/1_3_23_R1041_Duplex_HG002_1.fast5.tar"
    "s3://human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_Duplex/1_3_23_R1041_Duplex_HG002_2.fast5.tar"
    "s3://human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_Duplex/1_3_23_R1041_Duplex_HG002_3.fast5.tar"
)

for link in "${links[@]}"; do
	filename=$(basename "$link" | sed 's/.*\///')
	echo "unpacking ${filename%%.*}"
	aws s3 cp $link . --no-sign-request
	
	tar -xvf $filename
	mv 1_3_23_R1041_Duplex_HG002 "${filename%%.*}"
	rm -r $filename
done
