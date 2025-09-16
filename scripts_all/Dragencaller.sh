bwa_files="/home/kousis/work/ddrad/BWA_files_second"
NUM_CORES=28
REFERENCE="/home/kousis/work/ddrad/reference/T2T_ref.fasta"
output_vcf="/home/kousis/work/ddrad/vcf_second"
# Create chromosome list
output_file="chromosome_list.list"
grep "^>" "$REFERENCE" | awk '{print substr($1, 2)}' | sort | uniq | grep '^chr' > "$output_file"

pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/sorted_*.bam > tmp.txt
popd > /dev/null

mkdir -p tmp_v
mkdir -p $output_vcf


echo "Preprocessing Reference..."
#samtools faidx $REFERENCE
java -jar /home/kousis/work/tools/picard/build/libs/picard.jar  CreateSequenceDictionary \
       R=$REFERENCE
echo "Preprocessing Done"
Data_Pre_Processing() {
    local reference="$1"
    gatk ComposeSTRTableFile \
         -R "$reference" \
         -O "${output_vcf}/str_table.tsv"
    if [ $? -ne 0 ]; then
        echo "Error during Data pre processing"
        exit 1
    fi
}

Haplotype_caller() {
    local reference="$1"
    local input_bam="$2"
    local base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    samtools index "$input_bam"
    gatk CalibrateDragstrModel \
         -R "$reference" \
         -I "$input_bam" \
         -str "${output_vcf}/str_table.tsv" \
         -O "${output_vcf}/${base_name}_dragstr_model.txt"
    if [ $? -ne 0 ]; then
        echo "Error during Calibration"
        exit 1
    fi

    gatk HaplotypeCaller \
         -R "$reference" \
         -I "$input_bam" \
         -O "${output_vcf}/${base_name}.g.vcf" \
         -ERC GVCF \
         --tmp-dir tmp_v \
         --create-output-variant-index \
         --dragen-mode true \
         --dragstr-params-path "${output_vcf}/${base_name}_dragstr_model.txt"
    if [ $? -ne 0 ]; then
        echo "Error during call"
        exit 1
    fi

    # Uncomment if per-sample filtration is needed
    # gatk VariantFiltration \
    #      -V "${output_vcf}/${base_name}.g.vcf" \
    #      --filter-expression "QUAL < 50" \
    #      --filter-name "DRAGENHardQUAL" \
    #      --filter-expression "DP < 60" \
    #      --filter-name "LowDP" \
    #      -O "${output_vcf}/filtered_${base_name}.g.vcf" \
    #      --create-output-variant-index \
    #      --tmp-dir tmp_v
}

Merge() {
    local reference="$1"
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
         -L chromosome_list.list \
         --genomicsdb-workspace-path "${output_vcf}/my_database" \
         --tmp-dir tmp_v \
         --sample-name-map "${output_vcf}/cohort.sample_map" \
         --create-output-variant-index
    if [ $? -ne 0 ]; then
        echo "Error during GDB"
        exit 1
    fi
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
         -R "$reference" \
         -V gendb://${output_vcf}/my_database \
         -O "${output_vcf}/Merged.vcf" \
         --create-output-variant-index
    if [ $? -ne 0 ]; then
        echo "Error during call"
        exit 1
    fi

    # Select SNPs from the merged VCF
    gatk SelectVariants \
         -R "$reference" \
         -V "${output_vcf}/Merged.vcf" \
         --select-type-to-include SNP \
         -O "${output_vcf}/SNPs.vcf"
    if [ $? -ne 0 ]; then
        echo "Error during Variant Selection SNPs"
        exit 1
    fi

    # Select INDELs from the merged VCF
    gatk SelectVariants \
         -R "$reference" \
         -V "${output_vcf}/Merged.vcf" \
         --select-type-to-include INDEL \
         -O "${output_vcf}/INDELs.vcf"
    if [ $? -ne 0 ]; then
        echo "Error during Variant Selection Indels"
        exit 1
    fi
}

# Create cohort.sample_map file
create_sample_map() {
    while read -r bam_file; do
        base_name=$(basename "$bam_file" | sed 's/^sorted_//')
        echo -e "${base_name}\t${output_vcf}/${base_name}.g.vcf"
    done < "${bwa_files}/tmp.txt" > "${output_vcf}/cohort.sample_map"
}

# Data pre-processing
Data_Pre_Processing "$REFERENCE"

# Process each BAM file
count=0
while read -r bam_file; do
    running_jobs=$(jobs -p | wc -l)
    while [ "$running_jobs" -ge "$NUM_CORES" ]; do
        sleep 1
        running_jobs=$(jobs -p | wc -l)
    done
    Haplotype_caller "$REFERENCE" "$bam_file" &
    count=$((count + 1))
done < "${bwa_files}/tmp.txt"

wait

# Create sample map after all HaplotypeCaller jobs are done
create_sample_map

# Merge the GVCF files
Merge "$REFERENCE"

# Clean up temporary files
rm -r tmp_v
echo "All BAM files have been processed and merged."
