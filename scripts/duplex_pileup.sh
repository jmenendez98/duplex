#!/bin/bash

set -eux -o pipefail

duplex_bam=$1

/private/home/jmmenend/bin/modkit_v0.3.0 pileup-hemi \
	${duplex_bam} \
	-t 16 \
	--mod-threshold m:0.80 \
	--cpg \
	-r /private/groups/migalab/jmmenend/ref/hg002_v1.0/hg002v1.0.fasta \
	-o "${duplex_bam%.bam}.pileup.bed" \
	--only-tabs \
	--log "${duplex_bam%.bam}.pileup-hemi.log"
