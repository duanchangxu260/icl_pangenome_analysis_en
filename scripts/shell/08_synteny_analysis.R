#!/usr/bin/env Rscript
# =============================================================================
# ICL Gene Synteny Analysis
# =============================================================================
# Based on BLAST or MCScanX results, plot synteny around ICL gene loci.
#
# Usage:
#   # Method A: Use MCScanX output
#   Rscript 08_synteny_analysis.R --mcscanx output.collinearity
#
#   # Method B: Simple locus comparison diagram (BLAST + GFF only)
#   Rscript 08_synteny_analysis.R --blast best_hits.txt --gff-dir gffs/
#
# Output: PDF synteny diagram
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript 08_synteny_analysis.R [options]\n\n")
  cat("Options:\n")
  cat("  --gff-dir    DIR    GFF3 directory (for locus visualization)\n")
  cat("  --blast      FILE   BLAST result (-outfmt 6)\n")
  cat("  --mcscanx    FILE   MCScanX collinearity output\n")
  cat("  -o, --output FILE   Output PDF (default: icl_synteny.pdf)\n")
  cat("  --ref        STR    Reference species (default: Ath)\n")
  cat("  --flank      INT    Flanking gene count (default: 5)\n")
  quit("no")
}

GFF_DIR    <- parse_arg("--gff-dir", ".")
BLAST_FILE <- parse_arg("--blast", NULL)
MCSCANX    <- parse_arg("--mcscanx", NULL)
OUTPUT     <- parse_arg("-o", parse_arg("--output", "icl_synteny.pdf"))
REF_SPEC   <- parse_arg("--ref", "Ath")
FLANK      <- as.integer(parse_arg("--flank", "5"))

cat("\n========== ICL Synteny Analysis ==========\n")

# ----- Parse GFF, get gene coordinates -----
parse_gff_genes <- function(gff_dir) {
  gff_files <- list.files(gff_dir, pattern = "\\.gff3?$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
  if (length(gff_files) == 0) return(NULL)

  all_genes <- list()
  for (f in gff_files) {
    lines <- readLines(f, warn = FALSE)
    lines <- lines[!grepl("^#", lines)]
    if (length(lines) == 0) next

    dat <- read.table(text = lines, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "")
    colnames(dat) <- c("seqid", "source", "type", "start", "end",
                        "score", "strand", "phase", "attributes")

    gene_rows <- subset(dat, type == "gene")
    for (i in seq_len(nrow(gene_rows))) {
      attr <- gene_rows$attributes[i]
      gene_id <- regmatches(attr, regexpr("ID=([^;]+)", attr))
      gene_id <- sub("ID=", "", gene_id)
      gene_name <- regmatches(attr, regexpr("Name=([^;]+)", attr))
      gene_name <- sub("Name=", "", gene_name)

      all_genes[[length(all_genes) + 1]] <- data.frame(
        file     = basename(f),
        species  = gsub("(\\.gff3?|_ICL.*)$", "", basename(f)),
        seqid    = gene_rows$seqid[i],
        start    = gene_rows$start[i],
        end      = gene_rows$end[i],
        strand   = gene_rows$strand[i],
        gene_id  = ifelse(length(gene_id) > 0, gene_id, NA),
        gene_name= ifelse(length(gene_name) > 0, gene_name, NA),
        stringsAsFactors = FALSE
      )
    }
  }
  return(do.call(rbind, all_genes))
}

genes <- parse_gff_genes(GFF_DIR)

# ----- Method A: MCScanX results -----
if (!is.null(MCSCANX) && file.exists(MCSCANX)) {
  cat("Parsing MCScanX results:", MCSCANX, "\n")

  # MCScanX collinearity format: two genes per line
  col_data <- read.table(MCSCANX, stringsAsFactors = FALSE, fill = TRUE)
  cat("Collinear gene pairs:", nrow(col_data), "\n")

  # Simplified visualization: show collinearity blocks around ICL genes
  # Here we make a simple collinearity connection diagram

  # Filter blocks containing ICL
  icl_blocks <- col_data[grepl("ICL|icl|AT3G21720", apply(col_data, 1, paste, collapse=" ")), ]
  cat("Blocks containing ICL:", nrow(icl_blocks), "\n")

  pdf(OUTPUT, width = 10, height = 6)
  plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 10),
       xlab = "", ylab = "", axes = FALSE,
       main = "ICL Synteny Blocks (MCScanX)")

  if (nrow(icl_blocks) > 0) {
    n_blocks <- min(nrow(icl_blocks), 20)
    for (i in seq_len(n_blocks)) {
      y <- 9 - i * (8 / n_blocks)
      lines(c(1, 9), c(y, y), lwd = 2, col = scales::hue_pal()(n_blocks)[i])
      text(0.5, y, icl_blocks[i, 1], cex = 0.6, adj = 1)
      text(9.5, y, icl_blocks[i, 2], cex = 0.6, adj = 0)
    }
  } else {
    text(5, 5, "No synteny blocks containing ICL found\nPlease check MCScanX input data", cex = 1.2)
  }
  dev.off()

# ----- Method B: Locus comparison diagram -----
} else if (!is.null(genes) && nrow(genes) > 0) {
  cat("Plotting locus diagram using gene coordinates\n")

  # Find ICL genes
  icl_genes <- genes[grepl("ICL|icl|AT3G21720", genes$gene_name) |
                     grepl("ICL|icl|AT3G21720", genes$gene_id), ]

  if (nrow(icl_genes) == 0) {
    # If no explicitly labeled ICL genes found, use all genes
    icl_genes <- genes
    cat("Note: No genes explicitly labeled as ICL found, showing all genes\n")
  } else {
    cat("Found", nrow(icl_genes), "ICL gene(s)\n")
  }

  pdf(OUTPUT, width = 12, height = max(4, nrow(icl_genes) * 1.2))

  # Sort by species
  icl_genes <- icl_genes[order(icl_genes$species), ]
  icl_genes$y <- seq_len(nrow(icl_genes))

  p <- ggplot(icl_genes) +
    # Chromosome / scaffold
    geom_segment(aes(x = start - 5000, xend = end + 5000,
                     y = y, yend = y),
                 linewidth = 0.5, color = "grey70") +
    # Gene arrow
    geom_segment(aes(x = start, xend = end,
                     y = y, yend = y,
                     color = strand),
                 linewidth = 4, arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    # Gene label
    geom_text(aes(x = (start + end) / 2, y = y + 0.5,
                  label = paste(gene_name, gene_id, sep = "\n")),
              size = 2.5, vjust = 0, lineheight = 0.8) +
    scale_color_manual(values = c("+" = "steelblue", "-" = "coral"),
                       name = "Strand") +
    scale_y_continuous(breaks = icl_genes$y,
                       labels = paste(icl_genes$species, icl_genes$seqid, sep = "\n")) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 7, lineheight = 0.8),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Genomic coordinate (bp)", y = "",
         title = "ICL Gene Locus Comparison",
         subtitle = paste("Showing", nrow(icl_genes), "gene(s) | flanking +/-5kb"))

  print(p)
  dev.off()

} else {
  cat("Error: Neither MCScanX results nor GFF data available\n")
  cat("Please provide --gff-dir or --mcscanx\n")
  quit(status = 1)
}

cat("Synteny diagram saved:", OUTPUT, "\n")
cat("\n========== Done ==========\n")
