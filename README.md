# 🧬 EPIC v2 Methylation Array — Custom Annotation Package Builder

[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-green)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Builds a custom Bioconductor-style annotation package for the **Illumina
MethylationEPIC v2 (EPIC-8v2)** array directly from the Illumina-supplied
manifest CSV and a representative IDAT file. The output is a fully compatible
`IlluminaMethylationAnnotation` object that plugs into any standard `minfi`
workflow.

> **Role:** I was the primary data analyst responsible for designing and
> implementing this pipeline during my postdoctoral research on brain
> methylation profiling with EPICv2.

---

## Why This Pipeline Exists

The Illumina MethylationEPIC v2 array (EPIC-8v2-0_A2) is a newer platform than
the original 450K / EPIC v1 arrays. At the time of this analysis:

- No official Bioconductor annotation package existed for this exact manifest version
- The reference package (`IlluminaHumanMethylationEPICv2anno.20a1.hg38`) covers
  a different probe set than the physical chip used in our USC batch
- Physical bead addresses on our chips differed from the manifest, requiring
  chip-specific validation

This pipeline validates probe addresses against the actual IDAT files and
builds a chip-specific annotation object that correctly represents the probes
physically present on our arrays.

---

## Output Files

| File | Description |
|---|---|
| `Locations.rda` | Genomic coordinates: chr, pos (hg38), strand per probe |
| `Manifest.rda` | Probe design: addresses, sequences, Type I/II, dye channel |
| `Islands.UCSC.rda` | CpG island context: Island / Shore / Shelf / OpenSea |
| `SNPs.Illumina.rda` | Illumina-reported SNPs overlapping probe targets |
| `Other.rda` | All remaining annotation columns from the manifest |
| `SNPs.*.rda` | dbSNP overlap objects (one per input `*Single.rda`) |
| `<pkgName>.rda` | Combined `IlluminaMethylationAnnotation` object |
| `sessionInfo.txt` | R + package versions for reproducibility |

**Package name:** `IlluminaHumanMethylationEPICanno.ilm8A2.hg38`

---

## Pipeline Structure

```
00_setup.R              # packages, file paths, package name constants
01_manifest_loading.R   # parse manifest CSV; validate addresses vs IDAT
02_annotation_objects.R # build Locations, Manifest, Islands.UCSC, SNPs.Illumina, Other
03_snp_overlap.R        # intersect probes with dbSNP GRanges objects
04_build_package.R      # assemble IlluminaMethylationAnnotation; save all .rda files
```

Each script sources the previous one — run `04_build_package.R` to execute
the full pipeline end-to-end.

---

## Annotation Object Structure

```
IlluminaMethylationAnnotation
│
├── Locations       — chr / pos / strand  (GRanges-compatible DataFrame)
├── Manifest        — Name / AddressA / AddressB / ProbeSeqA / ProbeSeqB /
│                     Type / NextBase / Color
├── Islands.UCSC    — UCSC_RefGene_Group / UCSC_RefGene_Name /
│                     Islands_Name / Relation_to_Island
├── SNPs.Illumina   — SNP_ID / SNP_DISTANCE / SNP_MinorAlleleFrequency
├── Other           — remaining manifest columns (Methyl450_Enhancer, etc.)
└── SNPs.*          — dbSNP overlap results (one DataFrame per input file)
```

---

## Probe Type Reference

| Probe type | Design | Addresses | Colour channels |
|---|---|---|---|
| **Type I** | Two-colour, two-address | AddressA + AddressB | Red + Green |
| **Type II** | One-colour, one-address | AddressA only | Green only |
| **Control** | QC probes | Various | Various |
| **SNP I/II** | Genotyping (rs probes) | As above | As above |

Type I probes use two bead addresses: AddressA measures the methylated allele,
AddressB the unmethylated allele. A missing address means the corresponding
bead is absent from this chip and the probe must be dropped.

---

## Key Design Decisions

**Address validation against IDAT:**  
The manifest lists all probes designed for the platform, but individual chip
batches may be missing some bead types. We read the physical bead addresses
from a representative IDAT file and drop any probes whose AddressA or AddressB
is absent — preventing downstream NAs in methylation matrices.

**Column name harmonisation:**  
Illumina uses different column names (e.g. `Infinium_Design_Type`) from minfi's
expected names (e.g. `Type`). All renaming is centralised in `02_annotation_objects.R`.

**OpenSea recoding:**  
Probes with no CpG island relationship have an empty string in the raw manifest.
We recode these to `"OpenSea"` to match the minfi convention used by all
existing annotation packages.

**`SNP_DISTANCE` duplicate column:**  
The original code listed `SNP_DISTANCE` twice when building `SNPs.Illumina`.
This is silently corrected in `02_annotation_objects.R`.

---

## Usage in minfi Workflows

```r
# Load the custom annotation
load("IlluminaHumanMethylationEPICanno.ilm8A2.hg38.rda")

# Use with an RGChannelSet
RGSet <- read.metharray.exp(
  base       = "path/to/idats",
  annotation = c(array = "IlluminaHumanMethylationEPIC",
                 annotation = "ilm8A2.hg38")
)

# Extract annotation slots
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm8A2.hg38)
```

---

## Input Files Required

| File | Description |
|---|---|
| `EPIC-8v2-0_A2.csv` | Illumina manifest CSV (supplied with array kit) |
| `*_Grn.idat` | Any one Green channel IDAT from your batch |
| `objects/*Single.rda` | Pre-computed dbSNP GRanges (optional; for SNP overlaps) |

---

## Dependencies

```r
BiocManager::install(c(
  "minfi",
  "illuminaio",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "IlluminaHumanMethylationEPICmanifest",
  "GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics"
))
install.packages("dplyr")
```

---

## Citation

If you use this annotation in your research, please cite:

This Repo

---

## Author

**Hamed Abdollahi**
University of South Carolina
HA25@mailbox.sc.edu
---

## License

MIT License — see [LICENSE](LICENSE) for details.
