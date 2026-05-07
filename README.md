# phymapr 🌍🌳

**phymapr** is an R package designed to streamline phylogeographic analysis. It takes time-scaled phylogenetic trees (e.g., from BEAST or TreeTime) and genomic metadata, performs ancestral state reconstruction (ASR), and generates beautiful, animated geographic transmission maps.

## Features
* **End-to-End Pipeline:** Go from a `.nwk` file to an animated map in a single function call.
* **Dynamic Mapping:** Automatically calculates map boundaries with appropriate padding based on your geographic metadata.
* **Ancestral State Reconstruction:** Uses `fastAnc` to trace and estimate the geographic origins of ancestral lineages.
* **Dual Map Styles:** Choose between fast-rendering polygon global maps or sleek, dark-mode satellite tiles.

## Installation

Because `phymapr` relies on advanced phylogenetic visualization tools hosted on Bioconductor (specifically `ggtree` and `treeio`), you must install those dependencies using `BiocManager` before installing the package from GitHub.

Run the following code in your R console:

```R
# 1. Install BiocManager if you don't already have it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 2. Install the required Bioconductor dependencies
BiocManager::install(c("ggtree", "treeio"))

# 3. Install the 'remotes' package to download directly from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# 4. Install phymapr
remotes::install_github("ADHS-Taylor/phymapr")
```

## Quick Start

Here is a basic example of how to run the full transmission mapping pipeline. 

```R
library(phymapr)

# Run the master pipeline wrapper
results <- generate_phylo_transmission(
  timetree_file = "path/to/your_tree.tre", 
  metadata_file = "path/to/your_metadata.csv",
  col_specimen = "SPECIMEN_NUMBER", 
  col_date = "date",                
  col_location = "location",        
  map_style = "polygon" # Options: "polygon" (light/fast) or "tile" (dark/detailed)
)

# The function returns a list containing three objects:
# 1. tree_plot     (A ggtree ggplot object)
# 2. map_animation (A gganimate object)
# 3. mcc_data      (The cleaned, snapped coordinate data frame)

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
