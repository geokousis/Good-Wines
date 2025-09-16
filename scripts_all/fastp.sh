#!/bin/bash
output_dir="/home/kousis/work/ddrad/radtag/trimmed_new"
mkdir -p "$output_dir"
file_list="fix_to_trim_new.txt"

# Read the file list and process in pairs
while read -r r1_file && read -r r2_file; do
    base_name=$(basename "$r1_file" | sed 's/.1\.fq.*//')
    echo "in: $base_name"
    # Fastp
    echo "Processing: $r1_file and $r2_file"
    fastp -l 36 -i "$r1_file" -I "$r2_file" -o "${output_dir}/${base_name}_1_trimmed.fastq.gz" -O "${output_dir}/${base_name}_2_trimmed.fastq.gz"
    echo "Processed: $r1_file and $r2_file"
done < "$file_list"
