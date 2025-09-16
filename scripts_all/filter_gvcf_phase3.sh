#!/bin/bash

# Ensure correct usage
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bed> <input.vcf.gz>"
    exit 1
fi

BED_FILE=$1
VCF_FILE=$2
FILTERED_BED="filtered_regions.bed"
FILTERED_VCF="filtered_output.vcf.gz"
FILTERED_SNP_VCF="filtered_snps.vcf.gz"

# Check if files exist
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found!"
    exit 1
fi

if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file not found!"
    exit 1
fi

# Step 1: Filter BED file to keep only lines where the last column is 'True'
grep 'True$' "$BED_FILE" > temp_filtered.bed

# Step 2: Compute the 95th percentile threshold correctly
TOTAL_LINES=$(wc -l < temp_filtered.bed)
CUTOFF_LINE=$((TOTAL_LINES * 95 / 100))

PERCENTILE_VALUE=$(awk '{print $4}' temp_filtered.bed | sort -n | awk -v cutoff="$CUTOFF_LINE" 'NR == cutoff')

if [ -z "$PERCENTILE_VALUE" ]; then
    echo "Error: Could not determine 95th percentile value."
    exit 1
fi

# Step 3: Keep only rows where the fourth column is GREATER THAN OR EQUAL TO the 95th percentile
awk -v threshold="$PERCENTILE_VALUE" '$4 >= threshold' temp_filtered.bed > "$FILTERED_BED"
wc -l "$FILTERED_BED"
# Step 4: Filter the VCF file using bcftools (keep only variants in top 5% regions)
bcftools view -R "$FILTERED_BED" "$VCF_FILE" -Oz -o "$FILTERED_VCF"

# Step 5: Index the filtered VCF file
bcftools index "$FILTERED_VCF"

# Step 6: Extract only SNPs from the filtered VCF file
bcftools view -v snps "$FILTERED_VCF" -Oz -o "$FILTERED_SNP_VCF"

# Step 7: Index the filtered SNP VCF file
bcftools index "$FILTERED_SNP_VCF"

# Cleanup temporary file
rm temp_filtered.bed

echo "Filtering complete."
echo "Filtered VCF (Top 5% Regions): $FILTERED_VCF"
echo "Filtered SNP VCF (Only SNPs from Top 5%): $FILTERED_SNP_VCF"
