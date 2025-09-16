#!/bin/bash

start_time=$(date +"%Y-%m-%d %H:%M:%S")
figlet -c "BWA MEME"
echo -e "\nThis section of the script performs mapping (+indexing, sort, flagstat) on a given list of files using bwa-meme\n\n\n"

T_N=30
#reference="/home/kousis/work/ddrad/reference/all_chr_111113_V2.fasta"
reference="/home/kousis/work/ddrad/reference/T2T_ref.fasta"
ans="no"
file_list="tomap.txt"


total_lines=$(wc -l < "$file_list")

# BWA Indexing
if [ "$ans" == "yes" ]; then
  echo "Indexing started"
  bwa-meme index -a meme "$reference" -t "$T_N"
  ./build_rmis_dna.sh "$reference"
  
  if [ $? -ne 0 ]; then
    echo "Indexing failed"
    exit 1
  fi
  echo "Indexing done"
fi

# Initialize counter
counter=1
# BWA-MEME
while IFS= read -r line1 && IFS= read -r line2; do
    echo "Processing file $counter of $((total_lines / 2)): $line1 & $line2"

    # BWA MEME
    filename=$(basename "$line1" | sed 's/_[^_]*_trimmed\.fastq\.gz//')
    echo "BWA Mapping on $filename"
    bwa-meme mem -7 -M -t "$T_N" -R "@RG\tID:$filename\tSM:$filename\tPL:ILLUMINA" "$reference" "$line1" "$line2" > "${filename}.sam"
    if [ $? -ne 0 ]; then
      echo "BWA mapping failed for $filename"
      exit 1
    fi
    cowsay "BWA done"

    # SAM Files Processing
    samtools view -@ "$T_N" -b "${filename}.sam" > "${filename}.bam"
    samtools sort -@ "$T_N" "${filename}.bam" -o "sorted_${filename}.bam"
    samtools flagstat -@ "$T_N" "sorted_${filename}.bam" > "LOG_${filename}.txt"

    if [ $? -ne 0 ]; then
      echo "Samtools processing failed for $filename"
      exit 1
    fi

    # Remove temporary data
    rm "${filename}.sam" "${filename}.bam"

    counter=$((counter + 1))
done < "$file_list"

# File management
mkdir -p BWA_LOGS_second BWA_files_second
mv LOG_* BWA_LOGS_second
mv sorted_*.bam BWA_files_second

# Calculate and print the total execution time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))
echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds / 60) % 60)) minutes $((runtime_seconds % 60)) seconds"
