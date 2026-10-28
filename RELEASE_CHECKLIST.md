# v1.0.0 release checklist

## Local release preparation

- [x] Version chosen: `v1.0.0`; recorded in `VERSION`.
- [x] Public allowlist includes production code, docs, current figure, CI and synthetic fixtures.
- [x] Private configs/sheets, real data/results, Conda caches, draft figures and prompt logs stay local.
- [x] README, changelog, migration/method guides, release notes and colleague update prepared.
- [x] Final source/package checks and clean-package regression tests recorded in `STRUCTURAL_REVIEW.md`.
- [x] Local release commit, annotated tag, source archives and SHA-256 checksums prepared.

## Before remote publication

- [ ] Review the release diff and intended remote visibility.
- [ ] Select a distribution license with the repository owner/institution; no license has been selected automatically.
- [ ] Confirm maintainer/contact and citation details.
- [ ] Confirm hosted CI passes after pushing the release branch/tag.
- [ ] Complete representative biological validation, including a fresh-read run where applicable. The stopped cohort test is not a completed end-to-end run.
- [ ] Push the release commit/tag and publish the release notes/assets when authorized.

Preparation does not push or publish remotely. Existing local data and archived figures are retained.
