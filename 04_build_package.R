# =============================================================================
# 04_build_package.R
# Project: EPIC v2 Methylation Array — Custom Annotation Package Builder
# Author:  [Your Name]
# Description: Assemble all annotation objects into a single
#              IlluminaMethylationAnnotation S4 object and save everything
#              to disk in Bioconductor-standard .rda format.
#
#   Output files:
#     Locations.rda       — genomic coordinates per probe
#     Manifest.rda        — probe design (type, address, sequence)
#     SNPs.Illumina.rda   — Illumina-reported SNP overlaps
#     Islands.UCSC.rda    — CpG island and gene body relationships
#     Other.rda           — remaining annotation columns
#     SNPs.*.rda          — dbSNP overlap objects (one per *Single.rda input)
#     <pkgName>.rda       — the combined IlluminaMethylationAnnotation object
#                           (e.g. IlluminaHumanMethylationEPICanno.ilm8A2.hg38.rda)
#
#   The final .rda can be loaded into any minfi-based analysis with:
#     load("<pkgName>.rda")
#     RGSet <- read.metharray.exp(..., annotation = c(array = "IlluminaHumanMethylationEPIC",
#                                                     annotation = "ilm8A2.hg38"))
#
# Dependencies: Run 03_snp_overlap.R first.
# =============================================================================

source("03_snp_overlap.R")


# ── 1. Save Individual Annotation Objects ─────────────────────────────────────
# Each annotation object is saved as a separate .rda file with xz compression.
# This matches the structure of existing Bioconductor annotation packages
# (e.g. IlluminaHumanMethylationEPICv2anno.20a1.hg38) where each slot is a
# separately loadable file.

dir.create(PATHS$output_dir, showWarnings = FALSE, recursive = TRUE)

# Standard annotation objects
standard_objects <- c("Locations", "Manifest", "SNPs.Illumina",
                       "Islands.UCSC", "Other")

for (obj_name in standard_objects) {
  out_path <- file.path(PATHS$output_dir, paste0(obj_name, ".rda"))
  save(list   = obj_name,
       file    = out_path,
       compress = "xz")
  message("Saved: ", out_path)
}

# dbSNP overlap objects (variable number depending on inputs)
for (obj_name in snp.objects) {
  out_path <- file.path(PATHS$output_dir, paste0(obj_name, ".rda"))
  save(list   = obj_name,
       file    = out_path,
       compress = "xz")
  message("Saved: ", out_path)
}


# ── 2. Assemble IlluminaMethylationAnnotation Object ─────────────────────────
# IlluminaMethylationAnnotation is the top-level Bioconductor S4 class that
# registers all annotation slot names, the array type, and the genome build.
# When minfi loads an annotation package, it calls getAnnotation() on this
# object to retrieve the relevant slot (e.g. "Locations", "Islands.UCSC").
#
# objectNames: all slot names — the five standard plus any dbSNP overlaps
# annotation : named character vector identifying array + annotation + genome
# defaults   : slots loaded by default when getAnnotation() is called
# packageName: determines the R package namespace this object belongs to

all_anno_names <- c(standard_objects, snp.objects)

anno_obj <- IlluminaMethylationAnnotation(
  objectNames = all_anno_names,
  annotation  = ANNO_STR,
  defaults    = all_anno_names,   # all slots loaded by default
  packageName = PKG_NAME
)

# Assign to the expected package-level variable name
assign(PKG_NAME, anno_obj)

message("\nIlluminaMethylationAnnotation object created:")
print(anno_obj)


# ── 3. Save Combined Annotation Package Object ────────────────────────────────
# This single .rda is what users load to access the full annotation package
# without installing it from a Bioconductor repository.

pkg_rda_path <- file.path(PATHS$output_dir, paste0(PKG_NAME, ".rda"))
save(list    = PKG_NAME,
     file     = pkg_rda_path,
     compress = "xz")

message("\nPackage object saved: ", pkg_rda_path)


# ── 4. Session Information ─────────────────────────────────────────────────────
# Always record session info alongside annotation outputs to ensure
# reproducibility — annotation package versions affect downstream analysis.

session_path <- file.path(PATHS$output_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_path)
message("Session info saved: ", session_path)

message("\n=== Annotation package build complete ===")
message("Package name : ", PKG_NAME)
message("Output dir   : ", PATHS$output_dir)
message("Objects saved:")
for (f in list.files(PATHS$output_dir, pattern = "\\.rda$")) {
  size <- file.size(file.path(PATHS$output_dir, f))
  message(sprintf("  %-55s  %s KB", f, round(size / 1024)))
}
