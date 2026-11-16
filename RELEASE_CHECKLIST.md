# v1.0.1 release checklist

## Preparation and validation

- [x] Use a new patch version, `v1.0.1`, preserving the existing `v1.0.0` tag.
- [x] Include the reviewed HaplotypeCaller, duplicate-marking and coverage-display fixes and regression coverage.
- [x] Document Java job concurrency, CPU/memory budgets and optional per-job tuning in the README.
- [x] Preserve the public allowlist; private configurations, datasets, logs, caches and archived diagrams stay local.
- [x] Record current validation in `STRUCTURAL_REVIEW.md` and the release notes.
- [x] Retain the established maintainer and contact in the README.
- [x] Publication authorized by the maintainer on 16 November 2026.

## Publication workflow

A `VERSION` change pushed to `main` starts `.github/workflows/publish-release.yml`. It runs the reusable Snakemake CI workflow, then builds source archives from the checked-out commit and publishes the versioned GitHub release with a source manifest and SHA-256 checksums. Manual reruns are available through `workflow_dispatch`; an existing release is not overwritten.

Check the hosted workflow result and the release assets after publication. Source archives must correspond to the release tag and include the README and tests for that version.

## Distribution metadata

No distribution license has been selected or added. Citation details can be added by the maintainer later; they do not change the validation results recorded here.
