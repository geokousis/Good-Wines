#!/bin/bash

# Check for required arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <BED_FILE> <SNP_VCF.gz> <GVCF_VCF.gz>"
    exit 1
fi

# Input files
BED_FILE=$1
SNP_VCF=$2
GVCF_VCF=$3

# Output files
FILTERED_SNP_VCF="filtered_snp.vcf.gz"
FILTERED_GVCF_VCF="filtered_gvcf.vcf.gz"
SORTED_BED="sorted.bed"

# Step 1: Sort the BED file and ensure proper tab separation
awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' "$BED_FILE" | sort -k1,1 -k2,2n > "$SORTED_BED"

echo "Sorted BED file saved to $SORTED_BED"

# Step 2: Filter SNP VCF using sorted BED file
bcftools view -R "$SORTED_BED" "$SNP_VCF" -Oz -o "$FILTERED_SNP_VCF"
tabix -p vcf "$FILTERED_SNP_VCF"

echo "Filtered SNP VCF saved to $FILTERED_SNP_VCF"

# Step 4: Filter GVCF using BED file
bcftools view -R "$SORTED_BED" "$GVCF_VCF" -Oz -o "$FILTERED_GVCF_VCF"
tabix -p vcf "$FILTERED_GVCF_VCF"

echo "Filtered GVCF saved to $FILTERED_GVCF_VCF"
