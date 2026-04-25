# =============================================================================
# 01_manifest_loading.R
# Project: EPIC v2 Methylation Array — Custom Annotation Package Builder
# Author:  [Your Name]
# Description: Parse the Illumina-supplied EPIC v2 manifest CSV into the
#              minfi manifest structure, then validate probe addresses against
#              a representative IDAT file from the batch.
#
#   The manifest CSV (EPIC-8v2-0_A2.csv) contains:
#     - TypeI probes   : two-colour, two-address probes (AddressA + AddressB)
#     - TypeII probes  : one-colour, one-address probes (AddressA only)
#     - Control probes : quality control probes
#     - SNP probes     : rs-prefixed probes for genotyping / sample tracking
#
#   The IDAT validation step identifies probes whose physical addresses
#   are not present in the array scan — these must be dropped from the
#   annotation to avoid NAs downstream.
# =============================================================================

source("00_setup.R")


# ── 1. Parse Manifest CSV ─────────────────────────────────────────────────────
# minfi:::read.manifest.EPIC() reads the Illumina CSV and returns:
#   $manifest     — a data.frame of all probe metadata (one row per probe)
#   $manifestList — a named list of DataFrames split by probe type:
#                     TypeI, TypeII, TypeControl, TypeSnpI, TypeSnpII
# Note: ::: accesses minfi's internal (non-exported) function.

message("Reading EPIC v2 manifest CSV...")
maniTmp       <- minfi:::read.manifest.EPIC(PATHS$manifest)
anno          <- maniTmp$manifest       # full annotation data.frame
manifestList  <- maniTmp$manifestList   # probe-type split for manifest package

message("  Total probes in manifest: ", nrow(anno))
message("  Probe type breakdown:")
print(table(substr(anno$Name, 1, 2)))


# ── 2. Load IDAT Address List ─────────────────────────────────────────────────
# readIDAT() reads a single IDAT file and returns probe intensities plus
# metadata. We only need the MidBlock field, which lists the physical bead
# addresses present on this specific array chip.
# Any Grn.idat from the batch is sufficient — bead addresses are chip-level,
# not sample-level.

message("\nReading sample IDAT to extract physical bead addresses...")
epic         <- readIDAT(PATHS$idat_sample)
address_epic <- as.character(epic$MidBlock)
message("  Physical addresses on chip: ", length(address_epic))


# ── 3. Identify and Report Dropped Probes ─────────────────────────────────────
# Some manifest entries reference probe addresses not physically present on
# the array chip (e.g. probes added in later manifest revisions but absent
# from this batch's chips). These must be excluded.
#
# We check both AddressA (always present for all probe types) and AddressB
# (only present for Type I probes — empty string for Type II).
# A probe is dropped if ANY of its addresses is missing from the chip.

dropCpGs_B <- anno$Name[anno$AddressB != "" & !anno$AddressB %in% address_epic]
dropCpGs_A <- anno$Name[anno$AddressA != "" & !anno$AddressA %in% address_epic]
dropCpGs   <- unique(c(dropCpGs_A, dropCpGs_B))

message("\nProbes dropped (address not found on chip): ", length(dropCpGs))
if (length(dropCpGs) > 0) {
  message("  Prefix breakdown of dropped probes:")
  print(table(substr(dropCpGs, 1, 2)))
}


# ── 4. Build Manifest Package Object ─────────────────────────────────────────
# IlluminaMethylationManifest is the Bioconductor S4 class that stores the
# probe design information (addresses, sequences, types) used by minfi for
# preprocessing (background correction, dye normalisation, etc.).
# This object is analogous to IlluminaHumanMethylationEPICmanifest but for v2.

message("\nBuilding IlluminaMethylationManifest object...")
IlluminaHumanMethylationEPIC8v2manifest <- IlluminaMethylationManifest(
  TypeI        = manifestList$TypeI,
  TypeII       = manifestList$TypeII,
  TypeControl  = manifestList$TypeControl,
  TypeSnpI     = manifestList$TypeSnpI,
  TypeSnpII    = manifestList$TypeSnpII,
  annotation   = "IlluminaHumanMethylationEPIC"
)

message("  TypeI probes  : ",
        nrow(getProbeInfo(IlluminaHumanMethylationEPIC8v2manifest, type = "I")))
message("  TypeII probes : ",
        nrow(getProbeInfo(IlluminaHumanMethylationEPIC8v2manifest, type = "II")))
