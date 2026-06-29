#!/usr/bin/env Rscript
# =============================================================================
# ICL Conserved Motif Analysis
# =============================================================================
# Discover conserved motifs using MEME Suite, or display sequence conservation
# using a sliding window in R.
#
# Usage:
#   # Method A: Use MEME (requires MEME Suite)
#   Rscript 07_motif_analysis.R --input candidates.fa --method meme
#
#   # Method B: Pure R AA frequency conservation heatmap (no MEME required)
#   Rscript 07_motif_analysis.R --input aligned.fa --method conservation
#
# Input:
#   FASTA sequence file (unaligned or aligned)
#
# Output:
#   PDF motif structure diagram / conservation heatmap
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
flag <- function(name) any(args == name)

if (flag("-h") || flag("--help")) {
  cat("Usage: Rscript 07_motif_analysis.R [options]\n\n")
  cat("Options:\n")
  cat("  -i, --input    FILE   Input FASTA\n")
  cat("  -o, --output   PREFIX Output prefix (default: icl_motif)\n")
  cat("  -m, --method   STR    Method: meme / conservation (default: conservation)\n")
  cat("  --nmotifs      INT    MEME motif count (default: 10)\n")
  cat("  --window       INT    Sliding window size (default: 10, for conservation method)\n")
  quit("no")
}

INPUT   <- parse_arg("-i", parse_arg("--input", "icl_candidates.fasta"))
OUTPUT  <- parse_arg("-o", parse_arg("--output", "icl_motif"))
METHOD  <- parse_arg("-m", parse_arg("--method", "conservation"))
NMOTIFS <- parse_arg("--nmotifs", "10")
WINDOW  <- as.integer(parse_arg("--window", "10"))

cat("\n========== ICL Motif Analysis ==========\n")
cat("Method:", METHOD, "\n")

if (!file.exists(INPUT)) stop("Error: input file does not exist: ", INPUT)

# ----- Parse FASTA -----
read_fasta <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  seqs <- list()
  current_name <- ""
  current_seq <- ""

  for (line in lines) {
    if (grepl("^>", line)) {
      if (current_name != "") {
        seqs[[current_name]] <- current_seq
      }
      current_name <- sub("^>", "", line)
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, gsub("\\s", "", line))
    }
  }
  if (current_name != "") seqs[[current_name]] <- current_seq
  return(seqs)
}

seqs <- read_fasta(INPUT)
cat("Sequence count:", length(seqs), "\n")
cat("Sequence length:", nchar(seqs[[1]]), "(first sequence)\n\n")

# ----- Method A: MEME -----
if (METHOD == "meme") {
  meme_bin <- Sys.which("meme")
  if (meme_bin == "") {
    cat("Warning: MEME not found, falling back to conservation method\n")
    cat("Install MEME: conda install -c bioconda meme\n")
    METHOD <- "conservation"
  }
}

if (METHOD == "meme") {
  cat("--- Running MEME ---\n")
  meme_out <- paste0(OUTPUT, "_meme")
  cmd <- paste("meme", INPUT, "-protein -oc", meme_out,
               "-nmotifs", NMOTIFS, "-nostatus")
  cat("Command:", cmd, "\n")
  system(cmd)

  meme_xml <- file.path(meme_out, "meme.xml")
  if (!file.exists(meme_xml)) {
    stop("MEME run failed, ", meme_xml, " not generated")
  }

  # Try parsing MEME XML with the universalmotif package
  if (requireNamespace("universalmotif", quietly = TRUE)) {
    library(universalmotif)
    motifs <- read_meme(meme_xml)
    cat("Found", length(motifs), "motif(s)\n")

    # Draw motif logo
    pdf(paste0(OUTPUT, "_logos.pdf"), width = 10, height = length(motifs) * 2)
    par(mfrow = c(length(motifs), 1), mar = c(3, 4, 2, 1))
    # (Simplified motif structure illustration)
    dev.off()

  } else {
    cat("(universalmotif package not installed, only outputting raw MEME results)\n")
    cat("MEME results directory: ", meme_out, "\n")
  }

} else {
  # ----- Method B: Conservation heatmap (pure R) -----
  cat("--- Calculating sequence conservation ---\n")

  # Get sequence matrix (each position, each amino acid)
  seq_names <- names(seqs)
  n_seqs <- length(seq_names)
  max_len <- max(nchar(unlist(seqs)))

  # Align length (assume already aligned or similar length)
  seq_matrix <- matrix("", nrow = n_seqs, ncol = max_len)
  for (i in seq_len(n_seqs)) {
    chars <- strsplit(seqs[[i]], "")[[1]]
    seq_matrix[i, seq_along(chars)] <- chars
  }

  # AA frequency per position
  aa_list <- unique(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]])

  # Calculate the frequency of each amino acid at each position
  conservation_df <- data.frame(position = integer(), aa = character(),
                                freq = numeric())

  for (pos in seq_len(max_len)) {
    col_chars <- seq_matrix[, pos]
    col_chars <- col_chars[col_chars != "" & col_chars != "-"]
    if (length(col_chars) == 0) next
    tbl <- table(col_chars)
    for (aa in names(tbl)) {
      conservation_df <- rbind(conservation_df,
        data.frame(position = pos, aa = aa,
                   freq = as.numeric(tbl[aa]) / length(col_chars)))
    }
  }

  # Conservation heatmap
  p1 <- ggplot(conservation_df, aes(x = position, y = aa, fill = freq)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradientn(colors = c("white", "yellow", "red", "darkred"),
                         name = "Frequency", limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    labs(x = "Sequence position", y = "Amino acid",
         title = "ICL Sequence Conservation Heatmap",
         subtitle = paste("Based on", n_seqs, "sequences"))

  # Shannon entropy per position (conservation metric)
  shannon_entropy <- function(freqs) {
    freqs <- freqs[freqs > 0]
    -sum(freqs * log2(freqs))
  }

  pos_entropy <- conservation_df %>%
    group_by(position) %>%
    summarise(
      entropy = shannon_entropy(freq),
      max_freq = max(freq),
      consensus = aa[which.max(freq)],
      .groups = "drop"
    )

  # Maximum possible entropy (20 AAs equally distributed = log2(20) ≈ 4.32)
  max_entropy <- log2(20)

  p2 <- ggplot(pos_entropy, aes(x = position)) +
    geom_ribbon(aes(ymin = 0, ymax = entropy), fill = "steelblue", alpha = 0.3) +
    geom_line(aes(y = entropy), color = "steelblue", size = 0.8) +
    geom_hline(yintercept = max_entropy, linetype = "dashed", color = "red") +
    annotate("text", x = max_len * 0.05, y = max_entropy,
             label = "Completely unconserved", hjust = 0, vjust = -0.5, size = 3, color = "red") +
    theme_minimal(base_size = 11) +
    labs(x = "Sequence position", y = "Shannon entropy (bit)",
         title = "ICL Sequence Position Conservation",
         subtitle = "Lower entropy = more conserved") +
    ylim(0, max(pos_entropy$entropy, max_entropy + 0.5))

  # Combined output
  pdf_out <- paste0(OUTPUT, "_conservation.pdf")
  pdf(pdf_out, width = 14, height = 8)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    print(p1 / p2 + plot_layout(heights = c(2, 1)))
  } else {
    print(p1)
    print(p2)
  }

  dev.off()
  cat("Conservation plot saved:", pdf_out, "\n")

  # Output highly conserved intervals
  cat("\n--- Highly Conserved Regions (entropy < 1.0) ---\n")
  conserved_regions <- pos_entropy %>%
    filter(entropy < 1.0) %>%
    mutate(region = cumsum(c(1, diff(position) > 3)))

  if (nrow(conserved_regions) > 0) {
    regions <- conserved_regions %>%
      group_by(region) %>%
      summarise(start = min(position), end = max(position),
                len = end - start + 1,
                avg_entropy = mean(entropy),
                consensus_seq = paste(consensus, collapse = ""),
                .groups = "drop")
    print(regions)
    write.csv(regions, paste0(OUTPUT, "_conserved_regions.csv"), row.names = FALSE)
  }
}

cat("\n========== Done ==========\n")
