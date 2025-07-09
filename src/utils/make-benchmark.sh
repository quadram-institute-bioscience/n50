#!/bin/bash

set -euxo pipefail


# Configuration
OUTDIR="benchmark-data"
SIMSEQS_BIN="bin/n50_simseqs"

# Create output directory
echo "Creating benchmark data directory..."
mkdir -p ${OUTDIR}

# Check if n50_simseqs exists
if [ ! -f "${SIMSEQS_BIN}" ]; then
    echo "Error: ${SIMSEQS_BIN} not found!"
    echo "Please ensure n50_simseqs is built and available in bin/"
    exit 1
fi

echo "Generating benchmark datasets..."

# 1. Bacterial assembly (FASTA) - ~4Mb total genome
echo "  → Bacterial assembly (FASTA)..."
${SIMSEQS_BIN} --fasta -o ${OUTDIR} -p bacterial_ \
    1*3M 2*500K 5*100K 20*10K 100*1K 1000*100

# 2. SILVA-like ribosomal RNA database (FASTA) - ~150Mb total
echo "  → SILVA-like rRNA database (FASTA)..."
${SIMSEQS_BIN} --fasta -o ${OUTDIR} -p silva_ \
    50000*1200 30000*1500 20000*1800

# 3. Mammalian assembly (FASTA) - ~1.5Gb total genome
echo "  → Mammalian assembly (FASTA)..."
${SIMSEQS_BIN} --fasta -o ${OUTDIR} -p mammalian_ \
    23*50M 100*1M 1000*100K 10000*10K 50000*1K

# 4. Illumina short reads (FASTQ) - ~1.5Gb sequencing data
echo "  → Illumina short reads (FASTQ)..."
${SIMSEQS_BIN} --fastq -o ${OUTDIR} -p illumina_ \
    10000000*150

# 5. Nanopore long reads (FASTQ) - ~2Gb sequencing data  
echo "  → Nanopore long reads (FASTQ)..."
${SIMSEQS_BIN} --fastq -o ${OUTDIR} -p nanopore_ \
    100000*8K 50000*15K 20000*25K 5000*50K

# 6. Mixed short assembly contigs (FASTA)
echo "  → Mixed short contigs (FASTA)..."
${SIMSEQS_BIN} --fasta -o ${OUTDIR} -p contigs_ \
    1000*5K 5000*1K 10000*500 20000*200

# 7. PacBio HiFi reads (FASTQ) - ~1Gb high-quality long reads
echo "  → PacBio HiFi reads (FASTQ)..."
${SIMSEQS_BIN} --fastq -o ${OUTDIR} -p pacbio_ \
    500000*12K 300000*15K 100000*20K

# 8. Viral genomes collection (FASTA)
echo "  → Viral genomes (FASTA)..."
${SIMSEQS_BIN} --fasta -o ${OUTDIR} -p viral_ \
    1000*30K 5000*10K 10000*5K 20000*2K

echo ""
echo "Benchmark datasets generated successfully!"
echo "Output directory: ${OUTDIR}"
echo ""
echo "Dataset summary:"
echo "  - Bacterial assembly:    ~4 MB (1,126 sequences)"
echo "  - SILVA rRNA database:   ~150 MB (100,000 sequences)"  
echo "  - Mammalian assembly:    ~1.5 GB (61,123 sequences)"
echo "  - Illumina reads:        ~1.5 GB (10M reads)"
echo "  - Nanopore reads:        ~2 GB (175,000 reads)"
echo "  - Mixed contigs:         ~36 MB (36,000 sequences)"
echo "  - PacBio HiFi reads:     ~1 GB (900,000 reads)"
echo "  - Viral genomes:         ~400 MB (36,000 sequences)"
echo ""
echo "Total estimated size: ~6.5 GB"
