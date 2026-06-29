#!/bin/bash
# =============================================================================
# Batch extract CDS sequences from genome + GFF annotations (using gffread)
# =============================================================================
# Iterates over each species subdirectory in a pan-genome directory; if a CDS
# file is missing, it is automatically extracted with gffread.
#
# Prerequisites:
#   Each species directory must contain:
#     - Genome FASTA (*.fa / *.fasta / *.fna)
#     - Annotation file (*.gff / *.gff3 / *.gtf)
#
# Usage:
#   bash extract_cds_batch.sh <pan_genome_dir>
#   bash extract_cds_batch.sh /data1/duanchangxu/Pan-Genome
#
# Source: DeepSeek ICL Pan-Genome Analysis Pipeline
# =============================================================================

set -euo pipefail

WORK_DIR="${1:-.}"

if [ ! -d "$WORK_DIR" ]; then
    echo "Error: directory not found: $WORK_DIR"
    exit 1
fi

if ! command -v gffread &> /dev/null; then
    echo "Error: gffread is not installed"
    echo "Install: conda install -c bioconda gffread"
    echo "  or:  sudo apt install gffread"
    exit 1
fi

echo "== Batch extract CDS sequences =="
echo ""

extracted=0
skipped=0
failed=0

for genome_dir in "$WORK_DIR"/*/; do
    genome_name=$(basename "$genome_dir")

    # Check if CDS already exists
    existing=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.cds.fa" -o -name "*.cds.fasta" -o -name "*.cds.fna" \
        -o -name "*.cds.ffn" -o -name "*.cds" -o -name "cds.fasta" \) \
        | head -1)

    if [ -n "$existing" ]; then
        echo "⊙ $genome_name — already has CDS, skipping"
        skipped=$((skipped + 1))
        continue
    fi

    # Find genome and annotation files
    GENOME=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
        ! -name "*cds*" ! -name "*pep*" ! -name "*protein*" | head -1)
    GFF=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.gff" -o -name "*.gff3" -o -name "*.gtf" \) | head -1)

    if [ -z "$GENOME" ] || [ -z "$GFF" ]; then
        echo "✗ $genome_name — missing genome or annotation file, skipping"
        failed=$((failed + 1))
        continue
    fi

    # Check for CDS features in the GFF
    cds_features=$(grep -ci "[^#].*CDS" "$GFF" 2>/dev/null || echo 0)
    if [ "$cds_features" -eq 0 ]; then
        echo "✗ $genome_name — no CDS feature in GFF, trying -g extraction"
        # If a gene.fasta (cDNA) file exists, use it directly
        gene_fa=$(find "$genome_dir" -maxdepth 1 -type f -name "*.gene.fasta" | head -1)
        if [ -n "$gene_fa" ]; then
            echo "  → found gene.fasta, using directly as CDS"
            cp "$gene_fa" "${genome_dir}/cds.fasta"
            extracted=$((extracted + 1))
            continue
        fi
    fi

    # Extract CDS with gffread
    OUTPUT="${genome_dir}/cds.fasta"
    echo "→ $genome_name extracting..."

    if gffread -g "$GENOME" -x "$OUTPUT" "$GFF" 2>/dev/null; then
        count=$(grep -c "^>" "$OUTPUT" 2>/dev/null || echo 0)
        echo "  ✓ extraction complete: $count CDS sequences"
        extracted=$((extracted + 1))
    else
        # If gffread fails due to seqid mismatch, try falling back to gene.fasta
        gene_fa=$(find "$genome_dir" -maxdepth 1 -type f -name "*.gene.fasta" | head -1)
        if [ -n "$gene_fa" ] && [ ! -f "$OUTPUT" ]; then
            echo "  ⚠ gffread failed, falling back to gene.fasta as CDS"
            cp "$gene_fa" "$OUTPUT"
            count=$(grep -c "^>" "$OUTPUT" 2>/dev/null || echo 0)
            echo "  ✓ copy complete: $count CDS sequences"
            extracted=$((extracted + 1))
        else
            echo "  ✗ extraction failed"
            failed=$((failed + 1))
            rm -f "$OUTPUT"
        fi
    fi
done

echo ""
echo "== Done: $extracted newly extracted, $skipped skipped, $failed failed =="
echo ""
echo "Tip: extracted CDS files are named cds.fasta; consider merging them with merge_rename_cds.sh"

if [ "$failed" -gt 0 ]; then
    echo ""
    echo "Common failure reasons:"
    echo "  1. GFF seqid does not match FASTA header → troubleshoot with check_gff_fasta_match.sh"
    echo "  2. GFF is missing mRNA/CDS lines → check annotation completeness"
    echo "  3. Paths contain spaces/special characters → single quotes ' in directory names may cause issues"
fi
