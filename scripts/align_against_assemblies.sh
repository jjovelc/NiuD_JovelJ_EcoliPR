#!/usr/bin/bash

# Define paths
REFERENCE_GENOME="/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/NZ_HG941718.1.fasta"
WORKING_DIR="/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq"
BWA_INDEX="/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/NZ_HG941718.1"
TRIMMED_READS_DIR="/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq"
BAM_OUTPUT_DIR="/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/bamfiles"

# Validate input
R1_FILE=$1
if [[ -z "$R1_FILE" || ! -f "$R1_FILE" ]]; then
    echo "Error: Input R1 file not provided or does not exist."
    exit 1
fi

R2_FILE="${R1_FILE/_R1/_R2}"
if [[ ! -f "$R2_FILE" ]]; then
    echo "Error: Corresponding R2 file for $R1_FILE not found."
    exit 1
fi

# Extract sample name and prefix
SAMPLE_PREFIX=${R1_FILE%_R1.fastq.gz}
SAMPLE_NAME=$(basename "$SAMPLE_PREFIX")
SAMPLE_OUTPUT_FOLDER="$WORKING_DIR/$SAMPLE_NAME"

# Check if SAM file already exists
if [[ ! -f "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.sam" ]]; then
    echo "Processing sample: $SAMPLE_NAME"

    # Create output folders
    mkdir -p "$SAMPLE_OUTPUT_FOLDER" "$BAM_OUTPUT_DIR"

    # Read group string
    RG="@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:1_RG1_UNIT1\tLB:1_LIB1\tSM:${SAMPLE_NAME}"

    # Align reads
    echo "Aligning reads for sample $SAMPLE_NAME..."
    bwa mem -B 1 -M -t 8 -R "$RG" "$BWA_INDEX" "$R1_FILE" "$R2_FILE" > "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.sam" || { echo "Error: bwa mem failed"; exit 1; }

    # Convert SAM to BAM
    echo "Converting SAM to BAM..."
    samtools view -bh -F4 -T "$REFERENCE_GENOME" "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.sam" > "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.bam" || { echo "Error: samtools view failed"; exit 1; }

    # Sort and index BAM
    echo "Sorting and indexing BAM..."
    samtools sort -@ 8 "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.bam" -o "$BAM_OUTPUT_DIR/$SAMPLE_NAME.sorted.bam" || { echo "Error: samtools sort failed"; exit 1; }
    samtools index "$BAM_OUTPUT_DIR/$SAMPLE_NAME.sorted.bam" || { echo "Error: samtools index failed"; exit 1; }

    # Variant calling
    echo "Calling variants..."
    bcftools mpileup --threads 8 -Ou -f "$REFERENCE_GENOME" "$BAM_OUTPUT_DIR/$SAMPLE_NAME.sorted.bam" | \
        bcftools call --ploidy 1 -c -Ov -o "$BAM_OUTPUT_DIR/$SAMPLE_NAME.calls.vcf" || { echo "Error: bcftools call failed"; exit 1; }

    # Convert VCF to FASTQ
    echo "Converting VCF to FASTQ..."
    vcfutils.pl vcf2fq "$BAM_OUTPUT_DIR/$SAMPLE_NAME.calls.vcf" > "$BAM_OUTPUT_DIR/$SAMPLE_NAME.vcf2fq.fastq" || { echo "Error: vcfutils.pl failed"; exit 1; }

    # Convert FASTQ to FASTA
    echo "Converting FASTQ to FASTA..."
    seqtk seq -a "$BAM_OUTPUT_DIR/$SAMPLE_NAME.vcf2fq.fastq" > "$BAM_OUTPUT_DIR/$SAMPLE_NAME.fasta" || { echo "Error: seqtk failed"; exit 1; }

    # Clean up intermediate files
    echo "Cleaning up intermediate files..."
    rm "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.sam" "$SAMPLE_OUTPUT_FOLDER/$SAMPLE_NAME.bam"

else
    echo "Skipping sample: $SAMPLE_NAME (SAM file already exists)"
fi
