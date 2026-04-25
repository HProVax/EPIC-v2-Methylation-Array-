# =============================================================================
# 02_annotation_objects.R
# Project: EPIC v2 Methylation Array — Custom Annotation Package Builder
# Author:  Hamed Abdollahi
# Description: Harmonise the raw manifest columns to minfi's expected naming
#              conventions, then construct the five standard annotation objects
#              that together form the Bioconductor annotation package:
#
#   Locations    — genomic coordinates (chr, pos, strand) per probe
#   Manifest     — probe design (addresses, sequences, type, dye channel)
#   Islands.UCSC — CpG island / gene relationship annotations
#   SNPs.Illumina— Illumina-reported SNP overlaps near probes
#   Other        — all remaining annotation columns not in the above four
#
# Dependencies: Run 01_manifest_loading.R first.
# =============================================================================

source("01_manifest_loading.R")


# ── 1. Harmonise Column Names ─────────────────────────────────────────────────
# The Illumina CSV uses different column names from those expected by minfi.
# We map them here once so all downstream code uses minfi's conventions.
# The IlmnID column is removed because it will become the rowname.

anno$IlmnID <- NULL   # will be rownames, not a column

col_rename <- c(
  AddressA_ID           = "AddressA",
  AddressB_ID           = "AddressB",
  AlleleA_ProbeSeq      = "ProbeSeqA",
  AlleleB_ProbeSeq      = "ProbeSeqB",
  Infinium_Design_Type  = "Type",
  Next_Base             = "NextBase",
  Color_Channel         = "Color"
)

# Apply renaming: only rename columns that exist in the data
existing <- names(anno)
for (old_name in names(col_rename)) {
  new_name <- col_rename[[old_name]]
  if (old_name %in% existing) {
    existing[existing == old_name] <- new_name
  }
}
names(anno) <- existing

# Set rownames to the Illumina probe ID (IlmnID) for O(1) lookup
rownames(anno) <- anno$Name
message("Columns after harmonisation: ", paste(names(anno), collapse = ", "))


# ── 2. Filter to Locus Names Present in Existing EPIC Manifest ───────────────
# We restrict to probes whose CpG names (e.g. "cg00000029") are in the
# established IlluminaHumanMethylationEPICmanifest locus list.
# This ensures compatibility with existing minfi workflows that reference
# the canonical EPIC probe set.

locusNames <- getManifestInfo(IlluminaHumanMethylationEPICmanifest,
                              type = "locusNames")
anno_A     <- anno[anno$Name %in% locusNames, ]
message("Probes after filtering to EPIC locus names: ", nrow(anno_A))


# ── 3. Locations — Genomic Coordinates ───────────────────────────────────────
# Each probe maps to a single genomic position (the cytosine being measured).
# chr:  chromosome as "chr1", "chr2", …, "chrX", "chrY"
# pos:  integer genomic coordinate on hg38
# strand: "+" (forward) or "-" (reverse) — determines which strand the
#          cytosine is on; important for correctly interpreting methylation calls

Locations <- anno_A[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")

Locations$pos    <- as.integer(Locations$pos)
Locations$chr    <- paste0("chr", Locations$chr)   # add "chr" prefix for UCSC convention
Locations$strand <- ifelse(anno_A$Strand_FR == "F", "+", "-")

rownames(Locations) <- rownames(anno_A)
Locations <- as(Locations, "DataFrame")

message("Chromosome distribution:")
print(table(Locations$chr, exclude = NULL))


# ── 4. Manifest — Probe Design Information ────────────────────────────────────
# These columns are used by minfi for array preprocessing:
#   Name      — CpG probe identifier (e.g. "cg00000029")
#   AddressA  — physical bead address on the array (hexadecimal)
#   AddressB  — second address for Type I probes (empty for Type II)
#   ProbeSeqA — 50-mer probe sequence for allele A
#   ProbeSeqB — 50-mer probe sequence for allele B (Type I only)
#   Type      — "I" (two-colour) or "II" (one-colour)
#   NextBase  — base immediately 3' of the CpG (used for Type I extension)
#   Color     — dye channel: "Red" or "Grn" (Type I only)

Manifest <- anno_A[, c("Name", "AddressA", "AddressB",
                        "ProbeSeqA", "ProbeSeqB",
                        "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")
message("Manifest probe types: ", paste(names(table(Manifest$Type)), collapse = ", "))


# ── 5. Islands.UCSC — CpG Island and Gene Annotations ───────────────────────
# UCSC-derived CpG island and RefGene relationships:
#   UCSC_RefGene_Group   — gene body context (Body, 1stExon, TSS200, TSS1500, …)
#   UCSC_RefGene_Name    — gene symbol (e.g. "BRCA1")
#   UCSC_RefGene_Accession — RefSeq accession (e.g. "NM_007294")
#   Islands_Name         — CpG island identifier from UCSC
#   Relation_to_Island   — Island / N_Shore / S_Shore / N_Shelf / S_Shelf / OpenSea
#
# Note: probes with no island relationship have an empty string in the raw
# manifest; we recode these as "OpenSea" to match the standard minfi convention.

Islands.UCSC <- anno_A[, c("UCSC_RefGene_Group",
                            "UCSC_RefGene_Name",
                            "UCSC_RefGene_Accession",
                            "UCSC_CpG_Islands_Name",
                            "Relation_to_UCSC_CpG_Island")]

names(Islands.UCSC) <- c("UCSC_RefGene_Group",
                          "UCSC_RefGene_Name",
                          "UCSC_RefGene_Accession",
                          "Islands_Name",
                          "Relation_to_Island")

# Recode empty strings to "OpenSea" — standard minfi convention
Islands.UCSC$Relation_to_Island[Islands.UCSC$Relation_to_Island == ""] <- "OpenSea"
Islands.UCSC <- as(Islands.UCSC, "DataFrame")

message("CpG island relationship distribution:")
print(table(Islands.UCSC$Relation_to_Island, exclude = NULL))


# ── 6. SNPs.Illumina — Illumina-Reported SNP Overlaps ────────────────────────
# Illumina reports known SNPs (from dbSNP) that overlap each probe's target CpG.
# A SNP under the probe can confound methylation measurement — the probe may
# fail to hybridise correctly if the genotype differs from the reference.
#   SNP_ID                  — rs identifier of the overlapping SNP
#   SNP_DISTANCE            — distance from CpG to SNP (0 = at the CpG itself)
#   SNP_MinorAlleleFrequency— MAF; low-MAF SNPs may be less problematic

SNPs.Illumina <- anno_A[, c("SNP_ID",
                             "SNP_DISTANCE",
                             "SNP_MinorAlleleFrequency",
                             "SNP_DISTANCE")]   # SNP_DISTANCE listed twice in original
# Drop the inadvertent duplicate column
SNPs.Illumina <- SNPs.Illumina[, !duplicated(names(SNPs.Illumina))]
SNPs.Illumina <- as(SNPs.Illumina, "DataFrame")


# ── 7. Other — Remaining Annotation Columns ───────────────────────────────────
# All manifest columns not captured in the four objects above are stored in
# "Other" for completeness. Two renaming conventions are applied:
#   _NAME → _Name   (UCSC capitalisation convention)
#   X450k_Enhancer → Methyl450_Enhancer  (legacy array reference)

# Define all columns already used
used_cols <- c(
  names(Manifest),
  names(SNPs.Illumina),
  c("CHR", "MAPINFO", "Strand", "Chromosome_36", "Coordinate_36", "Genome_Build"),
  c("UCSC_RefGene_Group", "UCSC_RefGene_Name", "UCSC_RefGene_Accession",
    "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")
)

Other      <- anno_A[, setdiff(names(anno_A), used_cols)]
other_nams <- names(Other)
other_nams <- sub("_NAME", "_Name", other_nams)
other_nams[other_nams == "X450k_Enhancer"] <- "Methyl450_Enhancer"
names(Other) <- other_nams
Other <- as(Other, "DataFrame")

message("Other columns retained: ", ncol(Other))
message(paste(" ", names(Other), collapse = "\n"))
