#!/bin/bash
# =============================================================================
# Check which species in a pan-genome directory are missing CDS sequence files
# =============================================================================
# Usage:
#   bash check_cds.sh <pan_genome_dir>
#   bash check_cds.sh /data1/duanchangxu/Pan-Genome
# =============================================================================

set -euo pipefail

WORK_DIR="${1:-.}"

if [ ! -d "$WORK_DIR" ]; then
    echo "Error: directory not found: $WORK_DIR"
    exit 1
fi

echo "== Checking CDS file status =="
echo ""

missing=0
present=0

for genome_dir in "$WORK_DIR"/*/; do
    genome_name=$(basename "$genome_dir")

    # Find CDS files (supports multiple extensions)
    cds_file=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.cds.fa" -o -name "*.cds.fasta" -o -name "*.cds.fna" \
        -o -name "*.cds.ffn" -o -name "*.cds" -o -name "cds.fasta" \) \
        | head -1)

    if [ -n "$cds_file" ]; then
        count=$(grep -c "^>" "$cds_file" 2>/dev/null || echo "?")
        echo "✓ $genome_name — $(basename "$cds_file") ($count sequences)"
        present=$((present + 1))
    else
        # Check if genome + annotation are available for extraction
        has_genome=$(find "$genome_dir" -maxdepth 1 -type f \
            \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
            ! -name "*cds*" | head -1)
        has_gff=$(find "$genome_dir" -maxdepth 1 -type f \
            \( -name "*.gff" -o -name "*.gff3" -o -name "*.gtf" \) | head -1)

        if [ -n "$has_genome" ] && [ -n "$has_gff" ]; then
            echo "✗ $genome_name — no CDS (extractable: genome + annotation available)"
        else
            echo "✗ $genome_name — no CDS (missing genome or annotation, cannot extract)"
        fi
        missing=$((missing + 1))
    fi
done

echo ""
echo "== Summary: $present have CDS, $missing missing =="
