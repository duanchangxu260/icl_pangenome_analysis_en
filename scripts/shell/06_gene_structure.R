#!/usr/bin/env Rscript
# =============================================================================
# ICL Gene Structure Diagram (Exon-Intron)
# =============================================================================
# Extract ICL gene exon/intron structures from GFF3 files and plot a comparison
# diagram.
#
# Usage:
#   Rscript 06_gene_structure.R --gff-dir test_data --gene-ids genes.txt
#   Rscript 06_gene_structure.R -d gff_files/ -g icl_genes.txt -o icl_structure.pdf
#
# Input:
#   GFF3 files (one per species, placed in --gff-dir or current directory)
#   Gene ID list (one gene_id per line, corresponding to the Parent of mRNA in GFF)
#
# Output:
#   Gene structure comparison PDF
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript 06_gene_structure.R [options]\n\n")
  cat("Options:\n")
  cat("  -d, --gff-dir  DIR   GFF3 file directory (default: .)\n")
  cat("  -g, --gene-ids FILE  Gene ID list (one per line)\n")
  cat("  -o, --output    FILE  Output PDF (default: icl_gene_structure.pdf)\n")
  cat("  --width         NUM   PDF width (default: 10)\n")
  cat("  --height        NUM   PDF height (default: 6)\n")
  quit("no")
}

GFF_DIR  <- parse_arg("-d", parse_arg("--gff-dir", "."))
GENE_IDS <- parse_arg("-g", parse_arg("--gene-ids", NULL))
OUTPUT   <- parse_arg("-o", parse_arg("--output", "icl_gene_structure.pdf"))
FIG_W    <- as.numeric(parse_arg("--width", "10"))
FIG_H    <- as.numeric(parse_arg("--height", "6"))

cat("\n========== ICL Gene Structure Diagram ==========\n")

# ----- Step 1: Find all GFF3 files -----
gff_files <- list.files(GFF_DIR, pattern = "\\.gff3?$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)

if (length(gff_files) == 0) {
  stop("No GFF3 files found in ", GFF_DIR)
}
cat("Found", length(gff_files), "GFF3 file(s)\n")

# ----- Step 2: Parse GFF3, extract exon coordinates -----
parse_gff_exons <- function(gff_path) {
  if (!file.exists(gff_path)) return(NULL)

  # Skip comment lines when reading GFF
  lines <- readLines(gff_path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]

  if (length(lines) == 0) return(NULL)

  dat <- read.table(text = lines, sep = "\t", stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
  colnames(dat) <- c("seqid", "source", "type", "start", "end",
                      "score", "strand", "phase", "attributes")

  # Filter exon lines
  exons <- subset(dat, type == "exon")

  if (nrow(exons) == 0) return(NULL)

  # Extract gene name and transcript name
  extract_attr <- function(attr_str, key) {
    m <- regmatches(attr_str, regexpr(paste0(key, "=([^;]+)"), attr_str))
    if (length(m) == 0) return(NA)
    sub(paste0(key, "="), "", m)
  }

  exons$Parent  <- sapply(exons$attributes, extract_attr, "Parent")
  exons$gene_id <- sapply(exons$attributes, extract_attr, "ID")
  exons$gene_id <- ifelse(is.na(exons$gene_id),
                          sapply(exons$attributes, extract_attr, "gene_id"),
                          exons$gene_id)

  # Infer species name from file name
  exons$species <- gsub("(\\.gff3?|_ICL.*)$", "", basename(gff_path))
  exons$species <- gsub("[._]ICL$", "", exons$species)

  return(exons)
}

cat("Parsing GFF files...\n")
all_exons <- do.call(rbind, lapply(gff_files, function(f) {
  res <- parse_gff_exons(f)
  if (!is.null(res)) res$file <- basename(f)
  return(res)
}))

if (is.null(all_exons) || nrow(all_exons) == 0) {
  stop("Failed to extract exon information from any GFF")
}

cat("Extracted", nrow(all_exons), "exons in total\n")

# ----- Step 3: Normalize to relative coordinates -----
# For each gene, use the first exon's start as 0
all_exons <- all_exons %>%
  group_by(file, Parent) %>%
  mutate(
    gene_start = min(start),
    gene_end   = max(end),
    rel_start  = start - gene_start,
    rel_end    = end - gene_start,
    gene_len   = gene_end - gene_start
  ) %>%
  ungroup()

# ----- Step 4: Plot gene structure diagram -----
cat("Plotting...\n")

# Try gggenes; fall back to manual ggplot2 if not installed
has_gggenes <- requireNamespace("gggenes", quietly = TRUE)

if (has_gggenes) {
  library(gggenes)

  p <- ggplot(all_exons, aes(xmin = rel_start, xmax = rel_end,
                              y = file, forward = strand == "+")) +
    geom_gene_arrow(fill = "steelblue", color = "grey30", size = 0.5,
                    arrowhead_height = unit(4, "mm"),
                    arrowhead_width  = unit(3, "mm")) +
    geom_subgene_arrow(aes(xsubmin = rel_start, xsubmax = rel_end,
                           fill = "CDS"), color = NA) +
    facet_wrap(~ file, scales = "free", ncol = 1) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Relative position (bp)", y = "",
         title = "ICL Gene Structure (Exon-Intron)",
         subtitle = paste("Total", length(unique(all_exons$file)), "species"))

} else {
  # Fallback: manual ggplot2 drawing
  cat("(gggenes not installed, using manual ggplot2 drawing)\n")

  # Assign a unique y position to each gene
  gene_labels <- unique(all_exons$file)
  all_exons$y_pos <- match(all_exons$file, gene_labels)

  # Distinguish positive and negative strands
  all_exons$dir <- ifelse(all_exons$strand == "+", 1, -1)

  p <- ggplot(all_exons) +
    # Gene backbone line
    geom_segment(aes(x = 0, xend = gene_len, y = y_pos, yend = y_pos),
                 linewidth = 1.2, color = "grey40") +
    # Exon rectangle
    geom_rect(aes(xmin = rel_start, xmax = rel_end,
                  ymin = y_pos - 0.3, ymax = y_pos + 0.3),
              fill = "steelblue", color = "grey20", size = 0.3) +
    scale_y_continuous(breaks = seq_along(gene_labels),
                       labels = gene_labels) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Relative position (bp)", y = "",
         title = "ICL Gene Structure (Exon-Intron)",
         subtitle = paste("Blue rectangle = exon, grey line = intron | Total",
                          length(gene_labels), "gene(s)"))
}

ggsave(OUTPUT, p, width = FIG_W, height = FIG_H)
cat("Gene structure diagram saved:", OUTPUT, "\n")

# ----- Step 5: Output statistics table -----
cat("\n--- Gene Structure Statistics ---\n")
stats <- all_exons %>%
  group_by(file, Parent, species) %>%
  summarise(
    exon_count = n(),
    total_cds_len = sum(end - start + 1),
    gene_span = max(gene_len),
    strand = first(strand),
    .groups = "drop"
  )

print(as.data.frame(stats))

# Save statistics table
write.csv(stats, sub("\\.pdf$", "_stats.csv", OUTPUT), row.names = FALSE)
cat("\nStatistics table saved:", sub("\\.pdf$", "_stats.csv", OUTPUT), "\n")
cat("\n========== Done ==========\n")
