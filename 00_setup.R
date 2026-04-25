# =============================================================================
# 00_setup.R
# Project: EPIC v2 Methylation Array — Custom Annotation Package Builder
# Author:  [Your Name]
# Description: Load required Bioconductor and CRAN packages, and define all
#              file paths used across the annotation pipeline. Edit the PATHS
#              block below to match your local directory structure before running.
#
#   This pipeline builds a custom Bioconductor-style annotation package for the
#   Illumina MethylationEPIC v2 (EPICv2 / EPIC-8v2) array, equivalent to the
#   existing IlluminaHumanMethylationEPICv2anno.20a1.hg38 package but derived
#   directly from the Illumina-supplied manifest CSV and your own IDAT files.
# =============================================================================

# ── 1. Package Loading ────────────────────────────────────────────────────────

packages <- c(
  # Core methylation array tools
  "minfi",                                         # array QC, normalisation, annotation framework
  "illuminaio",                                    # low-level IDAT file reader
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",  # reference annotation (for comparison)
  "IlluminaHumanMethylationEPICmanifest",          # existing EPIC manifest (locus name reference)

  # Bioconductor infrastructure
  "GenomicRanges",   # GRanges: genomic coordinate representation
  "IRanges",         # IRanges: interval arithmetic
  "S4Vectors",       # DataFrame: Bioconductor data frames
  "BiocGenerics",    # getAnnotation(), common generics

  # Utilities
  "dplyr"
)

for (pkg in packages) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}


# ── 2. File Paths — edit these to match your system ───────────────────────────

PATHS <- list(

  # Illumina-supplied manifest CSV for EPIC v2 array
  manifest = "/media/hamed/HDD/Brain/Second_Round/idat_second_round/USC_Q3264_EPICv2/S2_BrainMeth_EPICv2/illumina/Annotation/EPIC-8v2-0_A2.csv",

  # One representative IDAT file from your batch (used to verify probe addresses)
  # Any Grn.idat from the batch will do — we only use the probe address list
  idat_sample = "/media/hamed/HDD/Brain/Second_Round/idat_second_round/USC_Q3264_EPICv2/USC_Q3264_EPIC_idats/208139170097_R01C01_Grn.idat",

  # Output directory for .rda annotation objects and final package .rda
  output_dir = "/media/hamed/HDD/Brain/Second_Round/idat_second_round/USC_Q3264_EPICv2/S2_BrainMeth_EPICv2",

  # Directory containing pre-computed SNP GRanges objects (*Single.rda files)
  snp_objects_dir = "objects"
)

# Array and annotation metadata — defines the output package name
ANNO_STR <- c(
  array        = "IlluminaHumanMethylationEPIC",
  annotation   = "ilm8A2",
  genomeBuild  = "hg38"
)

# Derived package name: IlluminaHumanMethylationEPICanno.ilm8A2.hg38
PKG_NAME <- sprintf("%sanno.%s.%s",
                    ANNO_STR["array"],
                    ANNO_STR["annotation"],
                    ANNO_STR["genomeBuild"])

message("Output package name: ", PKG_NAME)
message("Output directory   : ", PATHS$output_dir)
