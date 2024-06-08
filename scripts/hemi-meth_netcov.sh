#!/bin/bash

plus=$1
minus=$2
output=$3

# Extract necessary columns from both files
#cut -f 1,2,3,12 $1 > "${1%.bed}_cov.bedgraph"
#cut -f 1,2,3,12 $2 > "${2%.bed}_cov.bedgraph"

source /private/home/jmmenend/software/anaconda3/etc/profile.d/conda.sh

conda activate bedtools
# Merge the files based on positions using bedtools unionbedg
bedtools unionbedg -i $plus $minus > "${output}_combo_cov.bedgraph"
conda deactivate

# Perform the arithmetic operation using awk and save the result
awk 'BEGIN{OFS="\t"} 
{
    # If the columns from plus.bed or minus.bed are missing (.), set them to 0
    plus = ($4 == "." ? 0 : $4);
    minus = ($5 == "." ? 0 : $5);
    print $1, $2, $3, plus - minus
}' "${output}_combo_cov.bedgraph" > "${output}_net_cov.bedgraph"

echo "Processing complete. Results saved to "${output}_net_cov.bedgraph""
