#!/usr/bin/env python3
"""
analyze_missingness.py
This script reads a VCF file (compressed or uncompressed) using pysam,
counts missing genotype calls per sample, outputs the counts in a CSV file,
and produces a barplot with a table of the top 20 samples with highest missing counts.
A progress bar is used to track the counting process.
Usage: python3 analyze_missingness.py <vcf_file>
"""

import sys
import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

if len(sys.argv) < 2:
    print("Usage: {} <vcf_file>".format(sys.argv[0]))
    sys.exit(1)

vcf_file = sys.argv[1]

# Open VCF file using pysam
try:
    vcf = pysam.VariantFile(vcf_file)
except Exception as e:
    print(f"Error opening VCF file: {e}")
    sys.exit(1)

# Get sample names from VCF header
samples = list(vcf.header.samples)
missing_counts = {sample: 0 for sample in samples}

print("Processing VCF file to count missing genotypes...")

# Loop through each record with a progress bar.
# Note: The total number of records may be unknown, so tqdm will run in indeterminate mode.
for rec in tqdm(vcf.fetch(), desc="Counting records", unit="records"):
    for sample in samples:
        gt = rec.samples[sample].get("GT")
        # Count as missing if any allele in the genotype is None.
        if gt is None or any(allele is None for allele in gt):
            missing_counts[sample] += 1

# Create a DataFrame from the dictionary of missing counts
df = pd.DataFrame(list(missing_counts.items()), columns=["Sample", "MissingCount"])

# Save the missing counts to a CSV file.
df.to_csv("missing_counts.csv", index=False)
print("Missing counts saved in missing_counts.csv")

# Sort the DataFrame by MissingCount in descending order.
df_sorted = df.sort_values("MissingCount", ascending=False)

# Extract the top 20 samples.
top20 = df_sorted.head(20)
print("Top 20 samples with highest missing counts:")
print(top20.to_string(index=False))
