# ICL Pan-Genome Analysis

**Pan-genome isocitrate lyase (ICL) gene family identification and analysis pipeline for strawberry (*Fragaria*).**

This repository contains two tool sets:
- **ICL Screening** — HMMER → E-value filtering → sequence extraction → Arabidopsis BLAST → orthogroup screening
- **CDS Preparation** — CDS detection, batch extraction, merging and prefixing for pan-genome species

---

## Dependencies

```bash
# Required
sudo apt install ncbi-blast+

# Recommended (optional, awk fallback available)
conda install -c bioconda seqkit gffread
```

Python is not required.

---

## Directory Structure

```
icl_pangenome_analysis/
├── README.md
├── LICENSE
└── scripts/
    ├── check_cds.sh              # Check which species lack CDS
    ├── extract_cds_batch.sh      # Batch extract CDS (gffread)
    ├── merge_rename_cds.sh       # Merge CDS and add species prefix
    ├── 01_filter_hmm.sh          # HMMER E-value filtering
    ├── 02_extract_sequences.sh   # Protein sequence extraction
    ├── 03_blast_arabidopsis.sh   # Arabidopsis BLAST + filtering
    ├── 04_ogg_filter.R           # Orthogroup screening (R)
    ├── extract_fasta_by_ids.awk  # AWK FASTA extraction helper
    └── run_pipeline.sh           # One-click pipeline runner
```

---

## Workflow

### Prerequisite: CDS Preparation

```bash
# 1. Check which genomes lack CDS
bash scripts/check_cds.sh /data1/duanchangxu/Pan-Genome

# 2. Batch extraction (requires genome.fa + annotation.gff in each species directory)
bash scripts/extract_cds_batch.sh /data1/duanchangxu/Pan-Genome

# 3. Merge and add species prefix
bash scripts/merge_rename_cds.sh /data1/duanchangxu/Pan-Genome all_cds_merged.fasta
```

### ICL Gene Family Screening

**Input requirements:** Each species subdirectory must contain:
- `*.ICL.txt` — HMMER `--tblout` output
- `*.pep.fasta` or `*.fa` — proteome

```bash
# One-click run (four steps)
bash scripts/run_pipeline.sh /mnt/e/ogg arabidopsis_proteome.fasta

# Or run step by step
bash scripts/01_filter_hmm.sh /mnt/e/ogg 1e-5
bash scripts/02_extract_sequences.sh /mnt/e/ogg
bash scripts/03_blast_arabidopsis.sh /mnt/e/ogg ath.fa 1e-5 4 40 100
Rscript scripts/04_ogg_filter.R orthofinder --orthogroups Orthogroups.tsv --reference AT3G21720 -c best_hits.txt
```

### Step Descriptions

| Step | Script | Input | Output |
|------|------|------|------|
| E-value filtering | `01_filter_hmm.sh` | `*.ICL.txt` | `*_filtered_ids.txt` |
| Sequence extraction | `02_extract_sequences.sh` | ID list + protein FASTA | `*_ICL_candidates.fasta` |
| Arabidopsis BLAST | `03_blast_arabidopsis.sh` | candidate sequences + ath.fa | `best_hits.txt` |
| OGG screening | `04_ogg_filter.R` | best_hits + Orthogroups.tsv | `ogg_filtered_ids.txt` |

### Key Parameters

- HMMER E-value: `≤ 1e-5` (adjustable in commands)
- BLAST E-value: `≤ 1e-5`, identity `≥ 40%`, length `≥ 100aa`
- Arabidopsis ICL reference gene: `AT3G21720`

### Three OGG Screening Methods

1. **OrthoFinder** (recommended) — provide `Orthogroups.tsv`; the R script auto-matches
2. **eggNOG** — run `emapper.py` then use `04_ogg_filter.R` to match the target OGG
3. **RBH** — reciprocal best hits via two-way BLAST

---

## FAQ

**ID mismatch**: HMMER output IDs must match the first space-delimited word in FASTA headers. If inconsistent, use `seqkit grep -r` or `--match-by-prefix`.

**GFF/FASTA seqid mismatch**: First compare `grep "^>" genome.fa | head` and `cut -f1 annotation.gff | sort -u`, then use `sed` to fix the first column of the GFF.

**Paths with spaces/Chinese characters**: gffread is sensitive to special characters; use `sed 's/ /_/g'` to uniformly rename directories.

---

## References

- HMMER: Eddy SR (2011) *PLoS Comput Biol*
- BLAST+: Camacho C et al. (2009) *BMC Bioinformatics*
- OrthoFinder: Emms DM & Kelly S (2019) *Genome Biology*
- gffread: Pertea G & Pertea M (2020) *F1000Research*
