#!/usr/bin/env Rscript
###############################################################################
# run_methylkit.R
#
# Differential methylation analysis from Bismark cytosine reports using
# methylKit. Called automatically by methylation_pipeline.sh when
# RUN_DMR=true, or can be run standalone.
#
# Usage:
#   Rscript run_methylkit.R --cov_dir <dir> --group_a "s1 s2" \
#           --group_b "s3 s4" --context CpG --out_dir <dir>
###############################################################################

suppressMessages({
  if (!requireNamespace("methylKit", quietly = TRUE)) {
    stop("Package 'methylKit' not installed. Install via BiocManager::install('methylKit').")
  }
  library(methylKit)
  library(optparse)
})

# ---- CLI args ----
option_list <- list(
  make_option("--cov_dir", type = "character", help = "Dir with *.CX_report.txt / .cov files"),
  make_option("--group_a", type = "character", help = "Space-separated sample names, group A"),
  make_option("--group_b", type = "character", help = "Space-separated sample names, group B"),
  make_option("--context", type = "character", default = "CpG", help = "CpG, CHG, or CHH"),
  make_option("--out_dir", type = "character", default = "./dmr_results"),
  make_option("--min_coverage", type = "integer", default = 10),
  make_option("--qvalue", type = "double", default = 0.01),
  make_option("--diff", type = "double", default = 25, help = "Min methylation % difference")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

group_a <- strsplit(opt$group_a, "\\s+")[[1]]
group_b <- strsplit(opt$group_b, "\\s+")[[1]]
all_samples <- c(group_a, group_b)
treatment  <- c(rep(1, length(group_a)), rep(0, length(group_b)))

cat("Group A (treatment=1):", paste(group_a, collapse = ", "), "\n")
cat("Group B (treatment=0):", paste(group_b, collapse = ", "), "\n")

# Locate per-sample CX/cytosine report files produced by Bismark
find_report <- function(sample) {
  hits <- list.files(opt$cov_dir, pattern = paste0("^", sample, ".*CX_report.*txt(\\.gz)?$"),
                      full.names = TRUE)
  if (length(hits) == 0) stop("No cytosine report found for sample: ", sample)
  hits[1]
}
file_list <- lapply(all_samples, find_report)

# ---- Read in methylation calls ----
obj <- methRead(
  location   = as.list(file_list),
  sample.id  = as.list(all_samples),
  assembly   = "custom",
  treatment  = treatment,
  context    = opt$context,
  pipeline   = "bismarkCytosineReport",
  mincov     = opt$min_coverage
)

# ---- Filter + normalize by coverage ----
filtered <- filterByCoverage(obj, lo.count = opt$min_coverage, lo.perc = NULL,
                              hi.count = NULL, hi.perc = 99.9)
normalized <- normalizeCoverage(filtered)

# ---- Unite samples at common positions ----
meth <- unite(normalized, destrand = FALSE)
cat("Sites retained after uniting:", nrow(meth), "\n")

# ---- Sample correlation / clustering (QC) ----
pdf(file.path(opt$out_dir, "sample_correlation.pdf"))
getCorrelation(meth, plot = TRUE)
clusterSamples(meth, dist = "correlation", method = "ward.D2", plot = TRUE)
PCASamples(meth)
dev.off()

# ---- Differential methylation test ----
diff_meth <- calculateDiffMeth(meth, mc.cores = 1)

# All differentially methylated positions passing thresholds
dmp <- getMethylDiff(diff_meth, difference = opt$diff, qvalue = opt$qvalue)
write.csv(as.data.frame(dmp), file.path(opt$out_dir, "DMPs.csv"), row.names = FALSE)

hyper <- getMethylDiff(diff_meth, difference = opt$diff, qvalue = opt$qvalue, type = "hyper")
hypo  <- getMethylDiff(diff_meth, difference = opt$diff, qvalue = opt$qvalue, type = "hypo")
write.csv(as.data.frame(hyper), file.path(opt$out_dir, "DMPs_hyper.csv"), row.names = FALSE)
write.csv(as.data.frame(hypo),  file.path(opt$out_dir, "DMPs_hypo.csv"),  row.names = FALSE)

# ---- Tile into regions (e.g. 1kb windows) for DMRs ----
tiles <- tileMethylCounts(normalized, win.size = 1000, step.size = 1000, cov.bases = 3)
meth_tiles <- unite(tiles, destrand = FALSE)
diff_tiles <- calculateDiffMeth(meth_tiles, mc.cores = 1)
dmr <- getMethylDiff(diff_tiles, difference = opt$diff, qvalue = opt$qvalue)
write.csv(as.data.frame(dmr), file.path(opt$out_dir, "DMRs_1kb_tiles.csv"), row.names = FALSE)

cat("\nDone.\n")
cat(" DMPs (positions):", nrow(dmp), "\n")
cat(" Hyper-methylated:", nrow(hyper), " | Hypo-methylated:", nrow(hypo), "\n")
cat(" DMRs (1kb tiles):", nrow(dmr), "\n")
cat("Results written to:", opt$out_dir, "\n")
