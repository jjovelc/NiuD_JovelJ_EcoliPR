import glob, os, re  # Import libraries for file manipulation (glob), system commands (os), and regular expressions (re).

# Get a list of reference genome files (in .fasta format) from the specified directory.
referenceGenomes = glob.glob("/home/matthew/YanJan/RefGenomes/*.fasta")

# Iterate over each reference genome file.
for referenceGenome in referenceGenomes:
    # Define the output folder for this reference genome by replacing ".fasta" in its filename with an empty string.
    sampleOutputFolder = '/home/matthew/YanJan/RefAssemblies/' + re.sub('\.fasta', '', os.path.basename(referenceGenome))

    # Get a list of all paired-end read files ending with "R1.paired.fastq" from the specified directory.
    files = glob.glob('/home/matthew/YanJan/TrimmedReads/*R1.paired.fastq')

    # Define the folder for BAM files corresponding to this reference genome.
    bamOutputFolder = '/home/matthew/YanJan/BamFiles/' + re.sub('\.fasta', '', os.path.basename(referenceGenome))
    os.system('mkdir ' + bamOutputFolder)  # Create the BAM output folder if it doesn't exist.

    # Iterate over each paired-end read file.
    for file in files:
        # Define the sample prefix by removing "_R1.paired.fastq" from the file name.
        samplePrefix = re.sub('_R1\.paired\.fastq', '', file)

        # Extract the sample name from the sample prefix.
        sampleName = os.path.basename(samplePrefix)

        # Define the output folder for this sample.
        sampleOutputFolder1 = sampleOutputFolder + '/' + sampleName

        # Check if the SAM file for this sample already exists to avoid redundant processing.
        if not os.path.isfile(sampleOutputFolder1 + '/' + sampleName + '.sam'):
            print(sampleName)  # Print the sample name for logging.

            # Create the output folder for this sample.
            os.system('mkdir ' + sampleOutputFolder1)

            # Define file paths for paired-end reads (R1 and R2).
            sampleFastaR1 = samplePrefix + '_R1.paired.fastq'
            sampleFastaR2 = samplePrefix + '_R2.paired.fastq'

            print(f'=== Currently performing reference assembly for {sampleName} ===')  # Log the current sample.

            # Set the tab character for constructing the @RG string in BWA commands.
            tabC = "'\\'t"

            # Log the BWA command being executed.
            print(f'bwa mem -B 1 -M -t 8 -R @RG{tabC}ID:{sampleName}{tabC}PL:ILLUMINA{tabC}PU:1_RG1_UNIT1{tabC}LB:1_LIB1{tabC}SM:{sampleName} {referenceGenome} {sampleFastaR1} {sampleFastaR2} > {sampleOutputFolder1}/{sampleName}.sam')

            # Run BWA to align reads to the reference genome and output a SAM file.
            os.system(f'bwa mem -B 1 -M -t 8 -R @RG{tabC}ID:{sampleName}{tabC}PL:ILLUMINA{tabC}PU:1_RG1_UNIT1{tabC}LB:1_LIB1{tabC}SM:{sampleName} {referenceGenome} {sampleFastaR1} {sampleFastaR2} > {sampleOutputFolder1}/{sampleName}.sam')

            # Convert SAM to BAM, keeping only aligned reads.
            os.system(f'samtools view -bh -F4 -T {referenceGenome} {sampleOutputFolder1}/{sampleName}.sam > {sampleOutputFolder1}/{sampleName}.bam')

            # Sort the BAM file and save it to the BAM output folder.
            os.system(f'samtools sort {sampleOutputFolder1}/{sampleName}.bam > {bamOutputFolder}/{sampleName}.sorted.bam')

            # Index the sorted BAM file for downstream analysis.
            os.system(f'samtools index -b {bamOutputFolder}/{sampleName}.sorted.bam > {bamOutputFolder}/{sampleName}.sorted.bam.bai')

            # Perform variant calling using bcftools.
            os.system(f'bcftools mpileup -Ou -f {referenceGenome} {bamOutputFolder}/{sampleName}.sorted.bam | bcftools call --ploidy 1 -c -Ov -o {bamOutputFolder}/{sampleName}.calls.vcf')

            # Convert the VCF file to a FASTQ file using vcfutils.pl.
            os.system(f'vcfutils.pl vcf2fq {bamOutputFolder}/{sampleName}.calls.vcf > {bamOutputFolder}/{sampleName}.vcf2fq.fastq')

            # Convert the FASTQ file to a FASTA file using seqtk.
            os.system(f'seqtk seq -a {bamOutputFolder}/{sampleName}.vcf2fq.fastq > {bamOutputFolder}/{sampleName}.fasta')
