#!/bin/bash
# =============================================================================
# Merge all species CDS and add species prefix (prepare for downstream
# pan-genome analysis)
# =============================================================================
# Merge cds.fasta from each species subdirectory, adding a species name prefix
# to each sequence using seqkit. Fall back to awk when seqkit is unavailable.
#
# Usage:
#   bash merge_rename_cds.sh <pan_genome_dir> [output_file]
#   bash merge_rename_cds.sh /data1/duanchangxu/Pan-Genome all_cds.fasta
# =============================================================================

set -euo pipefail

WORK_DIR="${1:-.}"
OUTPUT="${2:-all_cds_merged.fasta}"

if [ ! -d "$WORK_DIR" ]; then
    echo "Error: directory not found: $WORK_DIR"
    exit 1
fi

echo "== Merging CDS sequences and adding species prefixes =="
echo ""

tmp_file=$(mktemp)
trap 'rm -f "$tmp_file"' EXIT

count=0
species_count=0

for genome_dir in "$WORK_DIR"/*/; do
    genome_name=$(basename "$genome_dir")

    # Find CDS file
    cds_file=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.cds.fa" -o -name "*.cds.fasta" -o -name "*.cds.fna" \
        -o -name "*.cds.ffn" -o -name "*.cds" -o -name "cds.fasta" \) \
        | head -1)

    if [ -z "$cds_file" ]; then
        echo "⊙ $genome_name — no CDS, skipping"
        continue
    fi

    # Use species name as prefix (replace spaces and special characters with underscores)
    prefix=$(echo "$genome_name" | sed 's/[^a-zA-Z0-9]/_/g' | sed 's/__*/_/g' | sed 's/_$//')

    seqs=$(grep -c "^>" "$cds_file" 2>/dev/null || echo 0)

    if command -v seqkit &> /dev/null; then
        seqkit rename -p "^($|.*)" -r "${prefix}_"'$1' "$cds_file" >> "$tmp_file"
    else
        # awk fallback: add prefix to > lines
        awk -v pfx="${prefix}_" '/^>/ { sub(/^>/, ">" pfx) } { print }' "$cds_file" >> "$tmp_file"
    fi

    echo "✓ $genome_name — $seqs sequences (prefix: ${prefix}_)"
    count=$((count + seqs))
    species_count=$((species_count + 1))
done

mv "$tmp_file" "$OUTPUT"

echo ""
echo "== Complete =="
echo "Species count: $species_count"
echo "Total CDS sequences: $count"
echo "Output file: $OUTPUT"
