# Data Preparation Guide

## Directory Structure

```
Pan-Genome/                        # Your pan-genome root directory
├── Fragaria_vesca_v1.0/
│   ├── Fvesca.ICL.txt             # HMMER hmmsearch --tblout output
│   ├── Fvesca.pep.fasta           # Proteome (or *.fa / *.protein.fasta)
│   ├── genome.fa                  # Genome assembly (for CDS extraction)
│   └── annotation.gff             # Annotation (for CDS extraction)
├── Fragaria_iinumae_v1.0/
│   └── ...
└── ...
```

## Key Notes

1. **HMMER output format**: Must be generated with `hmmsearch --tblout`. Default parsing reads column 1 (ID) and column 5 (E-value).
2. **FASTA header matching**: IDs from HMMER output must match the first space-delimited token in FASTA headers.
3. **Arabidopsis reference**: Download TAIR10 proteome; ICL reference gene is `AT3G21720`.
4. **Path naming**: Avoid spaces and non-ASCII characters. Recommend replacing with underscores.
