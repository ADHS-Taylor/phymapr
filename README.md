# phymapr 🌍🌳

**phymapr** is an R package for phylogeographic analysis and genomic epidemiology. It takes time-scaled phylogenetic trees (from BEAST, TreeTime, or Nextstrain) and sample metadata, performs ancestral state reconstruction, infers transmission pathways with confidence scoring, and generates animated and static geographic transmission maps.

## What's New in v0.2.0

- **Static cumulative map** — a single PNG showing all transmission pathways at once, much easier to interpret than the animated GIF alone
- **Improved import rendering** — distant/import arrows now render in black for clear visibility against the map background
- **Floating node fix** — internal nodes with unreliable ASR placements are filtered using a snap-distance threshold
- **Better pathway distinction** — thicker lines with solid/longdash/dotdash patterns and increased alpha
- **Date formatting** — GIF subtitle now shows readable dates (e.g., "Jan 06, 2026") instead of raw decimal years
- **Early pathway pruning** — `n_prune_early` removes deep evolutionary history from the animation

> To install the previous version: `remotes::install_github("ADHS-Taylor/phymapr@v0.1.0")`

## Installation

```R
# 1. Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Bioconductor dependencies
BiocManager::install(c("ggtree", "treeio"))

# 3. Install phymapr from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("ADHS-Taylor/phymapr")
```

## Quick Start

```R
library(phymapr)

results <- generate_phylo_transmission(
  timetree_file = "path/to/timetree.nexus",
  metadata_file = "path/to/metadata.tsv",
  col_specimen  = "accessionVersion",
  col_date      = "sampleCollectionDate",
  col_location  = "location",
  map_style     = "polygon",
  n_prune_early = 4
)

# Returns a list with 5 objects:
#   results$tree_plot     - ggtree phylogenetic tree (ggplot)
#   results$map_animation - animated transmission map (gganimate)
#   results$static_map    - cumulative static map (ggplot)
#   results$episodes      - local circulation episodes (data frame)
#   results$pathways      - transmission pathways with confidence (data frame)

# Save the static map (best for reports and presentations)
ggsave("transmission_map.png", results$static_map, width = 14, height = 10, dpi = 150)

# Save the animated GIF
library(gganimate)
anim <- animate(results$map_animation, width = 1200, height = 900, res = 120, fps = 10)
anim_save("transmission_map.gif", animation = anim)
```

## Example: Simulated Virus X, Y, and Z

The companion repository [phymap-workflow](https://github.com/ADHS-Taylor/phymap-workflow) includes three simulated datasets that demonstrate how data complexity affects transmission mapping.

### Virus X — Simple Domestic Spread

32 samples across 24 US cities over 10 months. A single introduction followed by clear sequential spread. The map shows clean, traceable direct pathways between locations.

| Tree | Transmission Map |
|------|-----------------|
| ![Virus X Tree](examples/virus_x_tree.png) | ![Virus X Map](examples/virus_x_map.gif) |

### Virus Y — Rapid Burst

30 samples from 22 cities in a compressed 2-month window. The short timeframe means many simultaneous events, producing a denser map where blended pacing separates overlapping transmission.

| Tree | Transmission Map |
|------|-----------------|
| ![Virus Y Tree](examples/virus_y_tree.png) | ![Virus Y Map](examples/virus_y_map.gif) |

### Virus Z — Clustered Transmission with Imports

33 samples concentrated in 12 cities over 2.5 months. Multiple introductions into the same locations create reticulate patterns with both direct pathways and import arrows, demonstrating the confidence classification system.

| Tree | Transmission Map |
|------|-----------------|
| ![Virus Z Tree](examples/virus_z_tree.png) | ![Virus Z Map](examples/virus_z_map.gif) |

## Key Features

- **Epidemiologic Inference Layer** — Collapses redundant local nodes into Local Residence Episodes and infers discrete transmission events between regions
- **Transmission Confidence Classification** — Direct Pathway (Likely), Indirect (Missing Intermediates), or Distant / Import based on SNP distance thresholds
- **Blended Non-Linear Pacing** — Smooths sampling gaps by blending calendar time with rank-order event progression
- **Auto-scaling Thresholds** — SNP distance cutoffs automatically calibrate to median branch divergence when data spans different evolutionary scales
- **Dynamic Geocoding** — Location strings resolved via OpenStreetMap with intelligent map bounds calculation
- **Multiple Output Formats** — Static PNG (cumulative), animated GIF, plus raw episode/pathway data frames

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `timetree_file` | — | Path to NEXUS/Newick time-scaled tree |
| `metadata_file` | — | Path to metadata CSV/TSV |
| `col_specimen` | `"SPECIMEN_NUMBER"` | Column matching tree tip labels |
| `col_date` | `"date"` | Collection date column |
| `col_location` | `"location"` | Geocodable location column |
| `map_style` | `"polygon"` | `"polygon"` (light) or `"tile"` (dark satellite) |
| `n_prune_early` | `4` | Remove N earliest pathways from animation |
| `threshold_direct_snp` | `2` | SNP cutoff for direct transmission |
| `threshold_indirect_snp` | `5` | SNP cutoff before import classification |
| `distant_rendering_style` | `"import_arrow"` | `"import_arrow"`, `"pulse"`, `"curve"`, or `"hide"` |
| `animation_pace_balance` | `0.5` | 1.0 = linear time, 0.0 = rank-based, 0.5 = blended |

## Full Pipeline (FASTA to Map)

For end-to-end use starting from raw aligned FASTA + metadata (tree building → molecular clock → transmission map), see the companion pipeline:

**[phymap-workflow](https://github.com/ADHS-Taylor/phymap-workflow)** — one command, no conda/Snakemake required:

```bash
python run_pipeline.py --fasta aligned.fasta.gz --metadata metadata.tsv.gz --prune
```

## Package Structure

| File | Role |
|------|------|
| `R/generate_phylo_transmission.R` | Master pipeline entry point |
| `R/data_prep.R` | Tree/metadata reading, geocoding, bounds |
| `R/tree_plotting.R` | Time-scaled tree visualization (ggtree) |
| `R/asr_modeling.R` | Brownian Motion ancestral state reconstruction |
| `R/epidemiologic_inference.R` | Episode collapsing, pathway scoring, pacing |
| `R/mapping.R` | Arc interpolation, static/animated map building |

## Citation

If you use phymapr in your work, please cite:

> Martins, T. (2026). phymapr: Phylogeographic Transmission Mapping in R. GitHub: https://github.com/ADHS-Taylor/phymapr

## License

MIT
