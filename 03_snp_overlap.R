# =============================================================================
# 03_snp_overlap.R
# Project: EPIC v2 Methylation Array — Custom Annotation Package Builder
# Author:  Hamed Abdollahi
# Description: Compute SNP overlaps with probe positions using GRanges.
#              For each pre-computed dbSNP GRanges object (*Single.rda), the
#              probe map is intersected to identify probes whose target CpG or
#              probe body overlaps a known SNP. This supplements the Illumina-
#              reported SNPs.Illumina object with population-level databases.
#
#   The probe position GRanges is built from Locations + Manifest, then
#   minfi:::.getProbePositionsDetailed() extends each probe to its full
#   hybridisation footprint (50 bp for the probe body).
#
#   Expected input: objects/*Single.rda files, each containing a GRanges
#   named grSnp<suffix> (e.g. grSnpCommon, grSnpAll).
#
# Dependencies: Run 02_annotation_objects.R first.
# =============================================================================

source("02_annotation_objects.R")


# ── 1. Build Probe Position GRanges ──────────────────────────────────────────
# Combine Locations (chr, pos, strand) with Manifest (Type) into a GRanges
# object. minfi:::.getProbePositionsDetailed() then extends each probe to its
# full 50 bp hybridisation footprint based on probe type and strand.
# This is necessary because a SNP anywhere in the probe body — not just at the
# CpG — can affect hybridisation efficiency and confound methylation readout.

message("Building probe position GRanges...")
map_df <- cbind(Locations, Manifest)

map <- GRanges(
  seqnames = map_df$chr,
  ranges   = IRanges(start = map_df$pos, width = 1),
  Strand   = map_df$strand,
  Type     = map_df$Type
)

# Extend probes to their full footprint
map <- minfi:::.getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

message("  Probe GRanges built: ", length(map), " probes")


# ── 2. Intersect with dbSNP GRanges Objects ───────────────────────────────────
# Each *Single.rda file in the objects/ directory contains a pre-computed
# GRanges of SNPs from a specific dbSNP build / frequency filter.
# Common file patterns:
#   grSnpCommon.Single.rda  — common SNPs (MAF > 0.01) from dbSNP
#   grSnpAll.Single.rda     — all SNPs regardless of frequency
#
# minfi:::.doSnpOverlap() performs the interval intersection and returns
# a DataFrame with overlap results for each probe.
# The resulting objects are named SNPs.<suffix> (e.g. SNPs.Common, SNPs.All).

snp.objects <- character(0)   # accumulates names of SNP overlap DataFrames

snp_files <- list.files(PATHS$snp_objects_dir, pattern = "*Single.rda",
                         full.names = FALSE)

if (length(snp_files) == 0) {
  warning("No *Single.rda files found in '", PATHS$snp_objects_dir, "'. ",
          "SNP overlap objects will be empty. ",
          "Place dbSNP GRanges .rda files in the objects/ directory to proceed.")
} else {
  message("Processing ", length(snp_files), " dbSNP GRanges file(s):")

  for (file in snp_files) {
    full_path       <- file.path(PATHS$snp_objects_dir, file)
    original_objname <- gsub("\\.rda$", "", file)     # e.g. "grSnpCommon"
    obj_name        <- gsub("grSnp", "SNPs.", original_objname)  # e.g. "SNPs.Common"

    message("  Loading: ", file, " → will create: ", obj_name)

    load(full_path)   # loads the grSnp* object into the environment

    # Perform probe–SNP overlap intersection
    overlap_result <- minfi:::.doSnpOverlap(map, get(original_objname))
    assign(obj_name, overlap_result)

    snp.objects <- c(snp.objects, obj_name)
    message("    Overlapping probes: ",
            sum(sapply(overlap_result, function(x) any(!is.na(x)))),
            " / ", length(map))
  }
}

message("\nSNP overlap objects created: ", paste(snp.objects, collapse = ", "))
