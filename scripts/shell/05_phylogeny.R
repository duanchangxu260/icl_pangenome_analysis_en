#!/usr/bin/env Rscript
# =============================================================================
# ICL gene family phylogenetic tree construction and visualization
# =============================================================================
# Input:  ICL candidate protein sequences FASTA (e.g., final_ICL_ogg.fa)
# Output:
#   - Multiple sequence alignment (aligned.fasta)
#   - Maximum likelihood tree (icl_tree.nwk)
#   - Phylogenetic tree figure PDF (icl_tree.pdf)
#
# Dependencies: mafft/muscle, iqtree/fasttree (command line)
#               R packages: ggtree, treeio, ape, ggplot2, ggmsa (optional)
#
# Usage:
#   Rscript 05_phylogeny.R --input final_ICL_ogg.fa --output icl_tree
#   Rscript 05_phylogeny.R -i candidates.fa -o icl_tree --aligner muscle
#   Rscript 05_phylogeny.R -i candidates.fa --skip-alignment (use existing alignment)
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(cowplot)
})

# ----- Argument parsing -----
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  # handle --name=value format
  for (a in args) {
    if (grepl(paste0("^", name, "="), a)) return(sub(paste0("^", name, "="), "", a))
  }
  return(default)
}

flag <- function(name) any(args == name)

if (flag("-h") || flag("--help")) {
  cat("Usage: Rscript 05_phylogeny.R [options]\n\n")
  cat("Options:\n")
  cat("  -i, --input FILE    Input FASTA (default: final_ICL_ogg.fa)\n")
  cat("  -o, --output PREFIX Output prefix (default: icl_tree)\n")
  cat("  --aligner CMD       Aligner: mafft (default) / muscle / clustalo\n")
  cat("  --tree-builder CMD  Tree builder: iqtree (default) / fasttree\n")
  cat("  --skip-alignment    Use existing alignment (input is already aligned FASTA)\n")
  cat("  --outgroup STR      Outgroup label (default: auto-detect Arabidopsis)\n")
  cat("  --threads INT       CPU threads (default: 4)\n")
  quit("no")
}

INPUT    <- parse_arg("-i", parse_arg("--input", "final_ICL_ogg.fa"))
OUTPUT   <- parse_arg("-o", parse_arg("--output", "icl_tree"))
ALIGNER  <- parse_arg("--aligner", "mafft")
TREEBLD  <- parse_arg("--tree-builder", "iqtree")
OUTGROUP <- parse_arg("--outgroup", NULL)
THREADS  <- parse_arg("--threads", "4")
SKIP_ALN <- flag("--skip-alignment")

ALIGNED  <- paste0(OUTPUT, "_aligned.fasta")
TREEFILE <- paste0(OUTPUT, ".nwk")
PDFFILE  <- paste0(OUTPUT, ".pdf")

cat("\n========== ICL Phylogenetic Tree Analysis ==========\n")
cat("Input sequences:", INPUT, "\n")
cat("Aligner:", ALIGNER, "\n")
cat("Tree builder:", TREEBLD, "\n")

if (!file.exists(INPUT)) {
  stop("Error: input file does not exist: ", INPUT)
}

# ----- Step 1: Multiple Sequence Alignment -----
if (!SKIP_ALN) {
  cat("\n--- Step 1: Multiple Sequence Alignment ---\n")

  aligner_cmd <- switch(ALIGNER,
    mafft   = paste("mafft --auto --thread", THREADS, INPUT, ">", ALIGNED),
    muscle  = paste("muscle -align", INPUT, "-output", ALIGNED),
    clustalo= paste("clustalo -i", INPUT, "-o", ALIGNED, "--threads", THREADS),
    stop("Unknown aligner: ", ALIGNER)
  )

  cat("Running:", aligner_cmd, "\n")
  system(aligner_cmd)

  if (!file.exists(ALIGNED) || file.info(ALIGNED)$size == 0) {
    stop("Alignment failed, no output file generated")
  }
  cat("Alignment complete:", ALIGNED, "\n")
} else {
  ALIGNED <- INPUT
  cat("Skipping alignment, using existing file:", ALIGNED, "\n")
}

# ----- Step 2: Build Phylogenetic Tree -----
cat("\n--- Step 2: Build Phylogenetic Tree ---\n")

if (TREEBLD == "iqtree") {
  tree_cmd <- paste("iqtree -s", ALIGNED, "-nt", THREADS, "-pre", OUTPUT, "-quiet")
  cat("Running:", tree_cmd, "\n")
  system(tree_cmd)
  # IQ-TREE outputs .treefile
  if (file.exists(paste0(OUTPUT, ".treefile"))) {
    file.copy(paste0(OUTPUT, ".treefile"), TREEFILE, overwrite = TRUE)
  }
} else if (TREEBLD == "fasttree") {
  tree_cmd <- paste("fasttree", ALIGNED, ">", TREEFILE)
  cat("Running:", tree_cmd, "\n")
  system(tree_cmd)
} else {
  stop("Unknown tree builder: ", TREEBLD)
}

if (!file.exists(TREEFILE) || file.info(TREEFILE)$size == 0) {
  stop("Tree building failed")
}
cat("Tree building complete:", TREEFILE, "\n")

# ----- Step 3: Read and Plot Phylogenetic Tree -----
cat("\n--- Step 3: Plot Phylogenetic Tree ---\n")

# Check if ggtree is installed
has_ggtree <- requireNamespace("ggtree", quietly = TRUE)

if (has_ggtree && requireNamespace("treeio", quietly = TRUE)) {
  # ----- Method A: ggtree (enhanced) -----
  library(ggtree)
  library(treeio)

  tree <- read.newick(TREEFILE)

  # Auto-detect outgroup (Arabidopsis)
  if (is.null(OUTGROUP)) {
    tip_labels <- tree$tip.label
    ath_tips <- grep("[Aa]th|AT[0-9]|Arabidopsis", tip_labels, value = TRUE)
    if (length(ath_tips) > 0) {
      OUTGROUP <- ath_tips[1]
    }
  }

  if (!is.null(OUTGROUP) && OUTGROUP %in% tree$tip.label) {
    tree <- root(tree, outgroup = OUTGROUP, resolve.root = TRUE)
    cat("Outgroup:", OUTGROUP, "\n")
  }

  # Main tree
  p <- ggtree(tree, ladderize = TRUE, size = 0.8) +
    geom_tiplab(size = 3, offset = 0.02, align = FALSE) +
    geom_tippoint(size = 2, color = "steelblue") +
    theme_tree2() +
    labs(title = "ICL Gene Family Phylogenetic Tree",
         subtitle = paste("Builder:", TREEBLD, "| Aligner:", ALIGNER))

  # Bootstrap values (if available)
  if (!is.null(tree$node.label) && any(tree$node.label != "")) {
    p <- p + geom_nodelab(aes(label = ifelse(as.numeric(label) >= 70, label, "")),
                          size = 2.5, vjust = -0.5)
  }

  ggsave(PDFFILE, p, width = 10, height = max(6, length(tree$tip.label) * 0.4))
  cat("Tree figure saved:", PDFFILE, "\n")

} else {
  # ----- Method B: ape base plot (no ggtree required) -----
  cat("(ggtree not installed, using ape base plot)\n")

  tree <- read.tree(TREEFILE)

  if (!is.null(OUTGROUP) && OUTGROUP %in% tree$tip.label) {
    tree <- root(tree, outgroup = OUTGROUP, resolve.root = TRUE)
  }

  pdf(PDFFILE, width = 10, height = max(6, length(tree$tip.label) * 0.4))
  plot(tree, cex = 0.8, no.margin = TRUE)
  title(main = "ICL Gene Family Phylogenetic Tree")
  add.scale.bar()
  # Add bootstrap values
  if (!is.null(tree$node.label) && any(tree$node.label != "")) {
    bs <- suppressWarnings(as.numeric(tree$node.label))
    bs[is.na(bs)] <- 0
    nodelabels(text = ifelse(bs >= 70, round(bs), ""),
               frame = "none", cex = 0.5, adj = c(1.2, -0.5))
  }
  dev.off()
  cat("Tree figure saved:", PDFFILE, "\n")
}

# ----- Step 4 (optional): MSA Visualization -----
if (has_ggtree && requireNamespace("ggmsa", quietly = TRUE)) {
  cat("\n--- Step 4: MSA Visualization ---\n")
  library(ggmsa)
  msa_pdf <- paste0(OUTPUT, "_msa.pdf")
  p_msa <- ggmsa(ALIGNED, start = 1, end = 200,
                 color = "Clustal", font = "DroidSansMono",
                 char_width = 0.5, seq_name = TRUE)
  ggsave(msa_pdf, p_msa, width = 14, height = max(4, length(tree$tip.label) * 0.3))
  cat("MSA figure saved:", msa_pdf, "\n")
}

cat("\n========== Done ==========\n")
cat("Output files:\n")
cat("  Alignment:", ALIGNED, "\n")
cat("  Tree file:", TREEFILE, "\n")
cat("  Tree figure:", PDFFILE, "\n")
