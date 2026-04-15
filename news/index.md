# Changelog

## ggNetView 0.1.0

### Initial CRAN release

- Initial submission of `ggNetView` to CRAN.
- Provides a unified, reproducible framework for analyzing and
  visualizing complex biological, ecological, and microbial association
  networks.
- Includes tools for building correlation and co-occurrence networks via
  `WGCNA`, `SpiecEasi`, `SparCC`, and standard correlation methods
  (Pearson, Spearman, Kendall).
- Computes node-level and network-level topological metrics, including
  robustness analyses.
- Supports module-level analyses (modularity detection, Zi-Pi
  classification, sample-level subgraph topology).
- Offers a large family of deterministic layout generators built on
  `ggraph` and `ggplot2` (bipartite, tripartite, quadripartite,
  pentapartite, circular-modules, petal, diamond, heart, star, and
  more), all with reproducible seeds.
- Ships 18 example datasets covering OTU tables, taxonomy tables,
  environmental metadata, PPI networks, and modularity examples.
