# phymapr 🌍🌳

**phymapr** is an R package designed to streamline phylogeographic analysis and genomic epidemiology. It takes time-scaled phylogenetic trees (e.g., from BEAST, TreeTime, or Nextstrain) and genomic metadata, performs ancestral state reconstruction (ASR), collapses local stays into circulation episodes, infers transmission pathways with multivariate confidence scoring, and generates beautiful, animated geographic transmission maps.

## Key Features
* **End-to-End Pipeline:** Go from a time-scaled tree file and metadata CSV to an animated map and phylogenetic tree in a single function call (`generate_phylo_transmission`).
* **Epidemiologic Inference Layer:** Collapses redundant local nodes into **Local Residence Episodes** and infers discrete transmission events between regions.
* **Transmission Confidence Classification:** Categorizes pathways into **Direct Pathway (Likely)**, **Indirect (Missing Intermediates)**, or **Distant / Import** based on SNP distance thresholds, temporal lag, and surveillance sampling density.
* **Blended Non-Linear Pacing:** Smooths out sampling gaps over time by blending true calendar time with rank-order event progression (`animation_pace_balance`).
* **Dynamic Mapping & Spatial Geocoding:** Automatically geocodes location strings (via OpenStreetMap) and calculates map bounding boxes with intelligent padding.
* **Dual Map Styles & Distant Rendering:** Choose between fast polygon world maps or dark satellite tile maps, with flexible import rendering (import arrows, pulses, bezier curves).

## Installation

Because `phymapr` relies on advanced phylogenetic visualization tools hosted on Bioconductor (specifically `ggtree` and `treeio`), install those dependencies using `BiocManager` before installing the package from GitHub.

Run the following code in your R console:

```R
# 1. Install BiocManager if you don't already have it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 2. Install required Bioconductor dependencies
BiocManager::install(c("ggtree", "treeio"))

# 3. Install the 'remotes' package
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# 4. Install phymapr from GitHub
remotes::install_github("ADHS-Taylor/phymapr")
```

## Quick Start

Here is a basic example of how to run the full transmission mapping pipeline:

```R
library(phymapr)

# Run the master pipeline wrapper
results <- generate_phylo_transmission(
  timetree_file = "path/to/your_tree.tre", 
  metadata_file = "path/to/your_metadata.csv",
  col_specimen = "SPECIMEN_NUMBER", 
  col_date = "date",                
  col_location = "location",        
  map_style = "polygon",                # Options: "polygon" or "tile"
  threshold_direct_snp = 2,            # Max SNP distance for direct transmission
  threshold_indirect_snp = 5,          # Max SNP distance before import classification
  distant_rendering_style = "import_arrow", # Options: "import_arrow", "pulse", "curve", "hide"
  animation_pace_balance = 0.5         # Blended pacing (1.0 = linear time, 0.0 = rank order)
)

# The function returns a list containing four objects:
# 1. tree_plot     (A ggtree ggplot object)
# 2. map_animation (A gganimate object)
# 3. episodes      (Data frame of local circulation residence episodes)
# 4. pathways      (Data frame of categorized transmission pathways)

# View the phylogenetic tree
print(results$tree_plot)

# Render and save the high-resolution animated map
library(gganimate)
high_res_map <- animate(
  results$map_animation, 
  width = 1200, 
  height = 900, 
  res = 120, 
  fps = 10
)

anim_save("transmission_map.gif", animation = high_res_map)
```

## Package Structure
- `R/data_prep.R`: Tree and metadata reading, ID normalization, geocoding, and spatial/temporal bounds calculation.
- `R/tree_plotting.R`: Time-scaled phylogenetic tree visualization using `ggtree`.
- `R/asr_modeling.R`: Brownian Motion Ancestral State Reconstruction (ASR).
- `R/epidemiologic_inference.R`: Local residence episode collapsing, transmission pathway confidence scoring, and non-linear blended pacing.
- `R/mapping.R`: Bezier arc, arrow, and pulse interpolation, ggplot map building, and gganimate transitions.
- `R/generate_phylo_transmission.R`: Master pipeline entry point.
