# v1.0.0 — Variant calling and genome analysis

Release candidate prepared on 28 October 2026.

This release brings optional population and local coverage analyses into the Snakemake calling workflow. A single config controls calling, filtering, PCA, phylogenetic trees, LEA structure analysis and coverage figures. Existing GVCF reuse and strict reference provenance remain in place.

## Main changes

- Separate PCA, TREE and STRUCTURE sections with independent enable switches; optional shared clade/rename/palette files.
- PCA/tree display lists select figure samples while preserving full-cohort analysis.
- LD pruning, NJ/UPGMA trees, sNMF K-range testing, cross-entropy summaries, Q matrices and admixture figures.
- Separate preparation/model/plotting jobs reuse unchanged work when display options or final K change.
- Optional per-BAM `CoveragePlots/` outputs using mosdepth fixed-window depth, optional GC correction and relative gain/loss screening.
- Coverage figures default to sequences of at least 50 kb, wrap after 16 sequences per row, and mark out-of-range values without removing analysis data.
- Public filtering handles missing annotations and scientific-notation INFO values; final-filter checks no longer depend on a pipe terminated by `head`.
- Updated workflow figure, documentation, tests and CI packaging.

## Upgrade

Use `config/config.yaml` as the full template. Existing flat configs remain supported; explicit legacy keys override nested equivalents. Do not mix duplicate settings in a new config. Relative paths resolve from the working directory.

Population and coverage analyses are disabled in the public template. Enable `postprocessing.enabled` and/or `coverage_plots.enabled` when needed. Activate the Snakemake environment before using `--use-conda`. Use `--rerun-incomplete` to resume interrupted work. Code/environment/filter changes can schedule downstream regeneration.

Private species-specific configs, biological data, results, environment caches and old diagram revisions are not part of the release; they remain local.

## Validation and limits

Local validation exercised 42 regression cases across two environments, plus clean-export structural/lint/DAG checks and real-tool synthetic integrations. See [the release review](../STRUCTURAL_REVIEW.md) for the exact results and environment-dependent skips. Tests use artificial data and installed tools; they do not replace evaluation of filters and coverage assumptions on representative biological samples. The full biological cohort run was stopped at the user's request, and is not reported as completed. Coverage regions are candidates, not absolute ploidy estimates or statistical CNV calls.

The candidate is prepared locally. Remote publication and selection of a distribution license remain maintainer actions; no license is invented by this update.
