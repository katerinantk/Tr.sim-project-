#!/usr/bin/env Rscript

# Compute omega scale factor for steric hindrance model
#
# Usage:
#   Rscript compute_omega_scale.R \
#     --rates    pause_escape_rates.csv \
#     --hkgenes  housekeeping_genes.txt \
#     --sample   SAMPLE_ID \
#     --out      scale_factor.csv
#
# Optionally provide --ref-median (the reference median omegaZeta for housekeeping
# genes from a matched published dataset) to also compute omega_scale_h.
# If not provided, only omega_scale_l is computed.
#
# The input rates CSV must have a chi column and either a gene_id or geneId column.
# It is expected to come from a pause-escape run (stericHindrance = FALSE).
#
# The housekeeping gene list should be a plain text file with one Ensembl gene ID
# per line (e.g. ENSG00000000419), or a CSV/TSV with a gene_id column.
#
# omega_scale_l = 0.2 / median(chi[housekeeping_genes])
# omega_scale_h = ref_median / median(chi[housekeeping_genes])
#
# Reference: Zhao, Liu & Siepel (2022) https://doi.org/10.1101/2022.10.19.512929

suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(
    prog = "compute_omega_scale.R",
    description = "Compute omega scale factor from pause-escape chi estimates"
)
parser <- add_option(parser, "--rates", type = "character", default = NULL,
    help = "Pause-escape rates CSV (must have chi column) [required]")
parser <- add_option(parser, "--hkgenes", type = "character", default = NULL,
    help = "Housekeeping gene list (one Ensembl ID per line, or CSV with gene_id column) [required]")
parser <- add_option(parser, "--sample", type = "character", default = NULL,
    help = "Sample ID to write in the output CSV [required]")
parser <- add_option(parser, "--out", type = "character", default = "scale_factor.csv",
    help = "Output CSV file [default: scale_factor.csv]")
parser <- add_option(parser, "--ref-median", type = "double", default = NULL,
    dest = "ref_median",
    help = "Reference median omegaZeta for housekeeping genes (for omega_scale_h) [optional]")
parser <- add_option(parser, "--append", action = "store_true", default = FALSE,
    help = "Append to existing output file instead of overwriting [default: FALSE]")

args <- parse_args(parser)

if (is.null(args$rates))   stop("--rates is required")
if (is.null(args$hkgenes)) stop("--hkgenes is required")
if (is.null(args$sample))  stop("--sample is required")

# ── Load rates ────────────────────────────────────────────────────────────────
rates <- read.csv(args$rates, stringsAsFactors = FALSE)

# Normalise column names: STADyUM uses geneId, UniMod uses gene_id
if ("geneId" %in% colnames(rates) && !"gene_id" %in% colnames(rates)) {
    colnames(rates)[colnames(rates) == "geneId"] <- "gene_id"
}

if (!"gene_id" %in% colnames(rates)) stop("rates file must have a gene_id or geneId column")
if (!"chi"     %in% colnames(rates)) stop("rates file must have a chi column")

cat(sprintf("Loaded %d genes from rates file\n", nrow(rates)))

# ── Load housekeeping genes ───────────────────────────────────────────────────
hk_raw <- readLines(args$hkgenes)
hk_raw <- hk_raw[nchar(trimws(hk_raw)) > 0 & !startsWith(trimws(hk_raw), "#")]

# If the first line looks like a header (no ENSG prefix), try reading as CSV
if (!any(grepl("^ENSG", hk_raw[1]))) {
    hk_tbl <- read.csv(args$hkgenes, stringsAsFactors = FALSE)
    gene_col <- intersect(c("gene_id", "geneId", "ensembl_gene_id"), colnames(hk_tbl))[1]
    if (is.na(gene_col)) stop("Cannot find gene ID column in housekeeping gene file")
    hk_genes <- hk_tbl[[gene_col]]
} else {
    hk_genes <- trimws(hk_raw)
}

cat(sprintf("Loaded %d housekeeping genes\n", length(hk_genes)))

# ── Intersect ─────────────────────────────────────────────────────────────────
hk_rates <- rates[rates$gene_id %in% hk_genes, ]
cat(sprintf("%d housekeeping genes found in rates table\n", nrow(hk_rates)))

if (nrow(hk_rates) < 10) {
    stop(sprintf(
        "Only %d housekeeping genes overlap with rates — too few for reliable calibration",
        nrow(hk_rates)
    ))
}

# ── Compute scale factors ─────────────────────────────────────────────────────
med_chi <- median(hk_rates$chi, na.rm = TRUE)
cat(sprintf("Median chi across housekeeping genes: %.6f\n", med_chi))

omega_scale_l <- 0.2 / med_chi
cat(sprintf("omega_scale_l (calibrated to 0.2 events/min): %.6f\n", omega_scale_l))

omega_scale_h <- NA_real_
if (!is.null(args$ref_median)) {
    omega_scale_h <- args$ref_median / med_chi
    cat(sprintf("omega_scale_h (calibrated to ref median %.4f): %.6f\n",
                args$ref_median, omega_scale_h))
} else {
    cat("omega_scale_h not computed (no --ref-median provided)\n")
}

# ── Write output ──────────────────────────────────────────────────────────────
result <- data.frame(
    sample_id      = args$sample,
    omega_scale_l  = omega_scale_l,
    omega_scale_h  = omega_scale_h,
    n_hk_genes     = nrow(hk_rates),
    median_chi_hk  = med_chi,
    stringsAsFactors = FALSE
)

if (args$append && file.exists(args$out)) {
    existing <- read.csv(args$out, stringsAsFactors = FALSE)
    # replace row if sample already exists, otherwise append
    existing <- existing[existing$sample_id != args$sample, ]
    result <- rbind(existing, result[, intersect(colnames(existing), colnames(result))])
}

write.csv(result, args$out, row.names = FALSE)
cat(sprintf("Written to %s\n", args$out))
