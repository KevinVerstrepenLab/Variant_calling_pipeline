# Short update for colleagues

**Subject: Variant-calling pipeline update — v1.0.0**

Hi everyone,

I have prepared an update to the Snakemake variant-calling pipeline. One configuration now controls the calling workflow and the optional downstream analyses:

- **PCA, phylogenetic trees and LEA structure analysis**, with shared clade information and separate sample lists for PCA/tree visualization.
- **Local coverage screening from sorted BAMs**, with per-sample figures in `CoveragePlots`, configurable window sizes, optional GC correction and candidate gain/loss regions.
- **Clearer coverage figures**, with robust y-axis limits, a default 50 kb display cutoff and up to 16 sequences per row.

The workflow continues to reuse validated existing GVCFs. Preparation, model fitting and plotting are separated, so changing figures does not unnecessarily repeat expensive analyses. Filtering was also hardened for missing annotations and scientific-notation values.

The README, workflow diagram, test suite and release packaging have been updated. Automated checks and synthetic integration tests cover the new branches; representative biological validation remains important before routine use. The optional analyses are disabled in the public template and can be enabled as needed.

Best,
