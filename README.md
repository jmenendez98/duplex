# Duplex Workflow:

### Data Pulled from:
```
https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_Duplex/1_3_23_R1041_Duplex_HG002_1.fast5.tar
https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_Duplex/1_3_23_R1041_Duplex_HG002_2.fast5.tar
https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_Duplex/1_3_23_R1041_Duplex_HG002_3.fast5.tar
```
### Basecalling: 
```
dorado-0.6.1 duplex sup,5mCG_5hmCG <pod5_folder> > <output.bam>
```
### Filtering:
Gives you a bam file only containing **paired** duplex reads
```
samtools view -b -h -d dx:1 <all_reads.bam> > <duplex.bam>
```
### Alignment:
```
methyl_bam=$1
output_prefix="$(basename ${methyl_bam%.bam})"

ref_fasta=$2
ref_name="$(basename -- $ref_fasta)"
ref_name="${ref_name%.*}"

temp="/private/groups/migalab/jmmenend/tmp/"

methyl_fastq="${output_prefix}.fastq"
samtools fastq -@ 64 -T MM,ML,Mm,Ml ${methyl_bam} > ${methyl_fastq}

conda activate /private/groups/migalab/jmmenend/.conda_envs/winnowmap
meryl count k=15 output merylDB_${ref_name}_k15 ${ref_fasta}
meryl print greater-than distinct=0.9998 merylDB_${ref_name}_k15 > ${ref_name}_repetitive_k15.txt
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
```
### Secphase:
```
NOT RUN YET
```

## Analyses:
### Duplex Rate Graphs:
Script, images and example log files from `dorado` are within [duplex_rate](duplex_rate/). The script uses log info from dorado to plot duplex rate by base on the left(as documented in the duplex portion of [dorado](https://github.com/nanoporetech/dorado)) and by read number on the right(which is contained in the log file).

### QC:
Nanoplot was run on all of the duplex files and zips for these can be found in [nanoplots](nanoplots/)
```
WIP
```

### Modkit:
Allows for aggregate analysis on hemi-methylation patterns.
```
WIP 
```

### 5mC bed tracks:
Allows us to look at hemi-methylation on the per-read basis
```
WIP
```