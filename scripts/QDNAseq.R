#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(QDNAseq)
  
  
})

# ---- 1. Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: QDNAseq_run.R <input_bam> <output_png>")
}

bamfile <- normalizePath(args[1])
out_png <- normalizePath(args[2])

# Derive useful paths
sample_name <- tools::file_path_sans_ext(basename(bamfile))
out_dir <- dirname(out_png)
seg_dir <- file.path(dirname(dirname(out_png)), "segmented_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(seg_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Run QDNAseq ----
message("Processing sample: ", sample_name)

bins <- getBinAnnotations(binSize = 1000, genome = "hg38")

readCounts <- binReadCounts(bins, bamfiles = bamfile, pairedEnds = TRUE)
readCountsFiltered <- estimateCorrection(readCounts)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers, method = "mean")
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)

# ---- 3. Plot ----
png(out_png, width = 1000, height = 600)
plot(copyNumbersSmooth)
dev.off()

# ---- 4. Export segments ----
exportBins(copyNumbersSegmented,
           file = file.path(seg_dir, paste0(sample_name, "_segments.txt")),
           type = "segments")

message("Done processing ", sample_name)
