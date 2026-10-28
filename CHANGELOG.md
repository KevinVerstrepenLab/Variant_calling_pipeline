# Changelog

## v1.0.0 — 2026-10-28

### Added

- A single sectioned run configuration, with compatibility for legacy flat configs.
- Optional PCA, NJ/UPGMA tree and LEA/sNMF structure analyses from the final SNP VCF.
- Independent PCA/tree display sample lists and shared sample/clade metadata.
- Separate preparation, model and plotting jobs with reusable internal checkpoints.
- Optional mosdepth coverage screening from new and existing sorted BAMs, with fixed windows, optional GC correction, relative gain/loss screening and PNG/PDF/QC exports.
- Coverage display controls: 50 kb default minimum sequence length, 16 sequences per row, and robust y limits with out-of-range markers.
- Updated workflow diagram, migration/method guides, release notes and a short colleague update.
- Public CI and synthetic regression fixtures, including real-tool filtering, BAM-to-plot and population-analysis checks.

### Fixed

- Guarded public hard-filter expressions against missing annotations and used decimal thresholds for scientific-notation INFO values.
- Replaced a pipefail-sensitive `bcftools | head` final-filter check with indexed record counts; forced safe replacement of intermediate VCF indexes on reruns.
- Preserved sNMF run numbers when cross-entropy contains nonfinite values.
- Prevented adoption of fitted models with missing or mismatching checkpoint provenance.
- Supported PCA figures without tree generation and tracked auto-selected best-K products explicitly.
- Matched LEA plot output names to the tested K range and included shared palette code in plotting signatures.
- Included CI and its artificial fixtures in the release allowlist while keeping private configurations, results and diagram history local.

### Retained

- Mixed paired/single-end processing, strict GVCF reference/sample provenance, indexed legacy-GVCF reuse, bounded interval imports and configurable SNP/INDEL filtering.

## Earlier development — 2026-10-12

- Consolidated the modular calling workflow and repaired repository-relative environment/script paths.
- Added mixed-read-layout support and explicit BAM/VCF index dependencies.
- Added SHA-256 reference/GVCF validation, persisted legacy adoption metadata and cohort manifests.
- Added directory-wide existing-GVCF reuse, uncompressed legacy inputs, bounded GenomicsDB batches and sample-map validation.
- Connected resource settings to tool threads and Java heaps, and added structural, provenance and DAG checks.
