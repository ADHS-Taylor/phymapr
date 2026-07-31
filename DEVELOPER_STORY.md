# phymapr 🌍🌳 Developer Story: Evolutionary Journey & Benchmarks

> [!NOTE]
> **Developer Journey Overview**  
> This report documents the complete evolutionary development of **phymapr**, tracing the pipeline from its initial working prototype, through the challenges encountered when porting to sampled real-world trees, the synthetic benchmarking on **Viruses X, Y, and Z**, and the final refactored **v0.2.0** pipeline release.

---

## 1. Stage 1: First Iteration (Original Prototype)

### 📖 The Initial Setup & Baseline
In the initial prototype phase, **phymapr** was designed as a lightweight set of R scripts to visualize time-scaled phylogenetic trees alongside continuous spatial ancestral state reconstruction (ASR). Using the original dataset and initial tree model, the pipeline successfully generated time-scaled `ggtree` plots and animated `gganimate` maps.

*Note: [Placeholder for user's initial work notes]*

### 🖼️ Artifacts: First Iteration Prototype

#### Time-Scaled Phylogenetic Tree
![Iteration 1 Tree Plot](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/iteration1_original_tree_plot.png)

#### Animated Geographic Transmission Map
![Iteration 1 Transmission Map GIF](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/iteration1_original_transmission_map.gif)

---

## 2. Stage 2: The Ported Version (Sampled Tree Issues)

### ⚠️ Cross-System Porting & Real-World Sampling Challenges
When porting the codebase to a secondary system and running it against real-world, sparsely sampled phylogenetic datasets (such as Pathoplexus genomic surveillance feeds), the visual output degraded significantly:

1. **Erratic Movement ("All Over the Place Too Quickly"):** With a sampled tree, continuous Brownian motion ASR produced erratic spatial trajectories. Lineages jumped rapidly across the map because linear frame animation forced long temporal gaps into rapid screen transitions.
2. **Deep Stem Distortion:** Unpruned early ancestor nodes created long ancestral stem branches that squished the recent, highly informative outbreak clades into a tiny fraction of the horizontal axis on the phylogenetic tree.
3. **Noisy Spatial Jitter:** Lineages remaining in a single location still generated micro-arcs across the map, creating visual clutter.

### 🖼️ Artifacts: Ported / Sampled Version (Unpruned & Unpaced)

#### Time-Scaled Tree with Long Early Stems (Unpruned)
![Ported Non-Pruned Tree Plot](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/iteration2_ported_non_pruned_tree_plot.png)

#### Animated Transmission Map with Rapid Sampling Jumps
![Ported Non-Pruned Transmission Map GIF](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/iteration2_ported_non_pruned_transmission_map.gif)

---

## 3. Stage 3: Synthetic Benchmark Exploration (Viruses X, Y, and Z)

To systematically diagnose, isolate, and solve these rate, sampling, and spatial diffusion challenges, three synthetic virus benchmark models were created: **Virus X**, **Virus Y**, and **Virus Z**.

---

### 🦠 Virus X: Linear Stepwise Wave Model

#### 🎯 Benchmark Objective
Simulate a clean, predictable geographic wave moving sequentially across the United States: East Coast $\rightarrow$ South $\rightarrow$ Midwest $\rightarrow$ West Coast $\rightarrow$ Pacific Northwest.

- **Dataset Size:** 32 samples across 24 US cities.
- **Mutational Profile:** Constant rate of 12 single-nucleotide substitutions per branch.
- **Key Insight:** Tested baseline spatial interpolation and ensured smooth temporal progression along a single geographical axis.

#### 🖼️ Artifacts: Virus X

##### Virus X Phylogenetic Tree
![Virus X Tree Plot](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/virus_x_tree_plot.png)

##### Virus X Animated Transmission Map
![Virus X Transmission Map GIF](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/virus_x_transmission_map.gif)

---

### 🦠 Virus Y: Hub-and-Spoke Airline & Ground Spread Model

#### 🎯 Benchmark Objective
Model dual-mode transmission: rapid long-distance jumps via major airport hubs (Atlanta, Chicago, LA, Dallas, Denver, Seattle, Miami) combined with localized ground diffusion.

- **Dataset Size:** 30 samples.
- **Mutational Profile:** 1–2 SNPs for rapid long-distance air travel branches; 4–7 SNPs for localized ground travel.
- **Key Insight:** Highlighted the need to visually distinguish long-distance import events from local regional circulation.

#### 🖼️ Artifacts: Virus Y

##### Virus Y Phylogenetic Tree
![Virus Y Tree Plot](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/virus_y_tree_plot.png)

##### Virus Y Animated Transmission Map
![Virus Y Transmission Map GIF](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/virus_y_transmission_map.gif)

---

### 🦠 Virus Z: Real-World Surveillance Imperfections Model

#### 🎯 Benchmark Objective
Stress-test the pipeline under severe real-world sampling anomalies and noisy data:
1. **Crowded Outbreak Clusters:** High-density sampling during local superspreading events (e.g. 9 near-identical genomes in Atlanta; 7 in LA).
2. **Unsampled Intermediates:** Large temporal and genetic gaps (e.g., 14 SNPs between Chicago and Denver with zero intervening samples).
3. **Unlinked Global Imports:** Highly drifted lineages (26 SNPs from root) appearing suddenly in Florida without domestic precursor nodes.
4. **Dead-End Singletons:** Lost transmission chains.

- **Key Insight:** Motivated the creation of the **Epidemiologic Inference Layer** (Local Residence Episodes + Transmission Confidence Categories) and **Blended Non-Linear Pacing**.

#### 🖼️ Artifacts: Virus Z

##### Virus Z Phylogenetic Tree
![Virus Z Tree Plot](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/virus_z_tree_plot.png)

##### Virus Z Animated Transmission Map
![Virus Z Transmission Map GIF](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/virus_z_transmission_map.gif)

---

## 4. Stage 4: The Final Redone Version (v0.2.0 Release)

### 🚀 Implementing the Solution Architecture
By integrating the insights gained from Viruses X, Y, and Z, the **phymapr v0.2.0** release was engineered to solve the sampled tree issues:

1. **Early Lineage Pruning (`n_prune_early = 4`):** Automatically trims long uninformative early stem branches, dramatically expanding the visual scale of recent transmission clades.
2. **Epidemiologic Inference Layer:** Snaps node locations to discrete spatial lookup entities and collapses regional stays into **Local Residence Episodes**, eliminating local spatial jitter.
3. **Multivariate Confidence Classification:**
   - 🟢 **Direct Pathway (Likely):** Low SNP distance ($\le 2$), short temporal window.
   - 🟡 **Indirect (Missing Intermediates):** Moderate SNP distance across sampling gaps.
   - 🔴 **Distant / Import:** Rendered clean import arrows/pulses for unlinked global imports.
4. **Blended Non-Linear Pacing (`animation_pace_balance = 0.5`):** Blends linear physical time with event rank-order weighting, smoothing out animation playback so sparse sampling periods no longer jump erratically.

### 🖼️ Artifacts: Final Redone Version (v0.2.0)

#### Pruned & Formatted Time-Scaled Tree
![Final Redone Tree Plot](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/final_redone_tree_plot.png)

#### Paced & Inferred Transmission Map Animation
![Final Redone Transmission Map GIF](file:///C:/Users/taylo/.gemini/antigravity-cli/brain/5a2303c2-0312-4689-be54-61292e807ccb/final_redone_transmission_map.gif)

---

## 5. Comparative Pipeline Evolution

| Feature / Metric | Stage 1: Prototype | Stage 2: Ported / Sampled | Stage 3: Viruses X, Y, Z | Stage 4: v0.2.0 Final Redone |
| :--- | :--- | :--- | :--- | :--- |
| **Package Structure** | Flat scripts | Flat scripts | Script benchmarks | Standard R package (`R/`, `NAMESPACE`) |
| **Stem Pruning** | None | None | None | `n_prune_early = 4` (Clades focused) |
| **Spatial Motion** | Continuous ASR jitter | Erratic jumps | Isolated in test models | Collapsed **Local Residence Episodes** |
| **Pacing Model** | Linear physical time | Linear (rapid gaps) | Stepwise / Hub / Noise | **Blended Non-Linear Pacing** ($\lambda = 0.5$) |
| **Pathway Categories** | Single line style | Single line style | Model-specific lines | **Direct**, **Indirect**, **Import Arrows** |
| **Visual Quality** | Baseline | Poor / Erratic | Benchmarked | **Production / High Resolution** |
