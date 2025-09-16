#!/bin/bash
# filter_vcf.sh
# This script takes:
#  1. A VCF file (can be compressed or uncompressed)
#  2. A text file (one sample per line) with samples to exclude.
#
# It then:
#   - Excludes the specified samples using bcftools view.
#   - Uses bcftools +fill-tags to add a missingness (MISS) tag,
#     and filters out sites with more than 5% missing genotypes.
#   - Extracts DP (depth) from the INFO field to compute the 5th percentile,
#     then filters out sites below that depth threshold.
#
# Usage: ./filter_vcf.sh <vcf_file> <exclude_samples.txt>

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <vcf_file> <exclude_samples.txt>"
    exit 1
fi

VCF=$1
EXCLUDE=$2

# Step 1: Exclude samples.
echo "Excluding samples listed in $EXCLUDE from $VCF..."
# The caret (^) before the file name tells bcftools to remove these samples.
bcftools view -S ^"$EXCLUDE" "$VCF" -Oz -o excluded_samples.vcf.gz --write-index

# Step 2: Filter sites with more than 10% missing genotypes.

echo "Filtering sites with missing rate > 10% (keeping sites with MISS <= 0.10)..."
bcftools view -i 'F_MISSING < 0.10' excluded_samples.vcf.gz -Oz -o filtered_missing.vcf.gz --write-index

# Step 3: Filter out sites with low depth.
echo "Extracting DP values from INFO field..."
bcftools query -f '%INFO/DP\n' filtered_missing.vcf.gz > dp_values.txt

# Compute the 5th percentile threshold of DP.
NUM=$(wc -l < dp_values.txt)
# Calculate the index (at least 1) corresponding to the 5th percentile.
PERCENTILE_INDEX=$(awk -v num="$NUM" 'BEGIN {printf "%d", (num*0.05 < 1 ? 1 : int(num*0.05)+1)}')
DP_THRESHOLD=$(sort -n dp_values.txt | sed -n "${PERCENTILE_INDEX}p")
echo "DP 5th percentile threshold: $DP_THRESHOLD"

echo "Filtering out sites with DP lower than threshold..."
bcftools filter -i "INFO/DP>=${DP_THRESHOLD}" filtered_missing.vcf.gz -Oz -o final_filtered.vcf.gz --write-index

echo "Filtering chr00..."
bcftools view final_filtered.vcf.gz | grep -v "chr00" > no_chr00.vcf
bgzip no_chr00.vcf
tabix -p vcf no_chr00.vcf.gz
echo "Splitting Multiallelic sites..."
bcftools norm -m -any no_chr00.vcf.gz  -o final_output.g.vcf.gz --write-index
echo "Final filtered VCF saved as final_output.g.vcf.gz"
