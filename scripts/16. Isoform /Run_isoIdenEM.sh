#!/bin/bash

# Define paths (replace with your actual paths)
fqinput="/path/to/DRS/01arrangeData"
extrac="/path/to/DRS/04Transcript/00extrac"
outpath="/path/to/DRS/04Transcript"

# Define reference files (replace with your actual files)
ref_path="/path/to/reference"
genome_fa="$ref_path/genome_rename.fa"
junc_bed="$ref_path/refdata_hg38.gtf2bed"
genes_gtf="$ref_path/genes/genes.gtf"

# Define software paths (replace with your actual paths)
seqkit_path="/path/to/seqkit"
minimap2_path="/path/to/minimap2"
samtools_path="/path/to/samtools"
transcriptclean_path="/path/to/TranscriptClean"
isoidenem_path="/path/to/IsoIdenEM.R"

# Samples to process
samples=("sample1" "sample2" "sample3")

# Create necessary directories
mkdir -p $outpath/01Align $outpath/02TranscriptClean $outpath/03reAlign $outpath/04FindIsoform

# Processing pipeline
for sample in "${samples[@]}"; do
    echo "Processing $sample..."

    # Extract reads from raw fastq
    $seqkit_path grep -f $extrac/${sample}_extra_readid.txt $fqinput/${sample}_DRS.fq.gz > $extrac/${sample}.fq

    # Align reads
    $minimap2_path -ax splice -t 20 --junc-bed $junc_bed $genome_fa $extrac/${sample}.fq | \
        $samtools_path sort -@ 20 | $samtools_path view -h - > $outpath/01Align/${sample}.sam

    # TranscriptClean run
    python $transcriptclean_path/TranscriptClean.py \
        --threads 20 \
        --sam $outpath/01Align/${sample}.sam \
        --genome $genome_fa \
        --spliceJns $extrac/${sample}_SpliceJns.txt \
        --variants $extrac/${sample}.vcf \
        --outprefix $outpath/02TranscriptClean/${sample} \
        --tmpDir $outpath/02TranscriptClean/

    # Re-align cleaned transcripts
    $minimap2_path -ax splice -t 20 --junc-bed $junc_bed $genome_fa $outpath/02TranscriptClean/${sample}_clean.fa | \
        $samtools_path sort -@ 20 | $samtools_path view -hb - > $outpath/03reAlign/${sample}_clean.bam

    # Run isoform identification
    Rscript $isoidenem_path \
        -c 30 \
        -g $genes_gtf \
        -r $genome_fa \
        -o $outpath/04FindIsoform/${sample}_ \
        --ignoreStrand TRUE \
        --MinReads2 5 --MinReads3 5 --MinReads4 5 --MinReads5 5 \
        $outpath/03reAlign/${sample}_clean.bam

    echo "$sample processing complete."
done

echo "All samples processed."
