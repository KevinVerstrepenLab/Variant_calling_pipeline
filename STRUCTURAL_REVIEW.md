# v1.0.0 pipeline and release review

Review date: 28 October 2026. Version: `1.0.0`. Scope: production rules, shared helpers, Python/R scripts, environment declarations, public configuration, test fixtures, documentation and the release allowlist.

## Findings corrected

| Finding | Resolution and evidence |
| --- | --- |
| CI and its required test files were excluded while diagram drafts/prompts were included | Public allowlist now includes CI and artificial fixtures; private configs, data/results, archives and diagram history remain local |
| Final filtering used `bcftools view | head` under `pipefail` | Replaced with indexed record counts; real final-filter rules passed for 5,000 SNP records and an empty INDEL VCF |
| Intermediate tabix indexes could survive a rerun and prevent replacement | Final-filter indexing uses `tabix -f` |
| Missing INFO annotations and integer thresholds could break/bypass compound hard filters | Public SNP/INDEL expressions guard each annotation and use floating-point thresholds; GATK tests passed with scientific notation and missing fields |
| Dropping nonfinite cross-entropy entries could renumber sNMF runs | Nonfinite values retain their run positions as NA; all-invalid K values stop explicitly |
| Public docs referenced a private run config and nonexistent distributed bundle | Examples now use the portable config and document that standalone/local archives are excluded |
| Structural checks referenced the wrong sample-template filename and missed nested rule directives | Corrected the filename, included release-required files/R source scanning, and inspected indented rules |

The new analysis branches were reviewed for explicit input/output dependencies, phase-specific parameters, metadata and sample-list handling, checkpoint invalidation and resource use. Existing GVCF validation, cohort maps, interval imports, genotyping, selection/filtering, concatenation and final filtering were also inspected. No full biological cohort job was restarted for this review.

## Validation performed

- Main local suite: **42 cases discovered, 40 passed and 2 skipped**. The two skipped bcftools genotype-masking cases passed separately in the Snakemake environment: all 42 cases were therefore exercised successfully across the two environments.
- Real-tool synthetic integration: GATK SNP/INDEL filters; empty and nonempty final-filter rules; samtools/mosdepth BAM-to-figure execution; PLINK PCA/tree preparation; R tree/PCA display subsets; LEA model and plotting execution.
- Reuse checks: changing K selection or sample lists retained matching models/preparation; seed changes refitted LEA; coverage display changes retained depth results; changed external BAM timestamps scheduled recomputation.
- Coverage regression checks: zero-depth windows, local gains/losses, GC trends, masked-window boundaries, display-size filtering, outlier retention in tables and multi-row figures.
- Snakemake 9.5.1 lint passed. The artificial mixed paired/single-end calling dry-run schedules 29 jobs.
- Dependency-free structural checks and staged whitespace checks passed. Python/R syntax, YAML, Markdown links and staged-file privacy are checked again during clean-package validation.
- Clean staged export: the same 42-case suite passed with the same two environment-dependent skips; both masking cases passed separately. Structural checks, lint and the 29-job dry-run also passed from that export, without private configs, standalone bundles or caches. Logs are retained with the local release artifacts.

Tool integrations use existing installed environments. Re-solving every Conda environment from scratch, hosted CI execution and a new complete FASTQ-to-all-analyses run are not claimed.

## Scope and remaining limitations

- The full biological cohort test was stopped at the user's request; it is not a completed end-to-end validation. All 187 coverage-plot jobs completed before stopping, while the calling cohort rebuild was still in progress.
- Calling filters require evaluation for each species/ploidy/library. Current final filtering retains samples and logs a warning if the sample-removal threshold would remove the entire cohort.
- A nonempty sample sheet/discovery result is currently required even for a fully reused GVCF cohort. Calling-stage paths should avoid whitespace.
- Display sample lists do not refit PCA/tree inference; full inferred Newick trees are retained. Unmatched clade IDs are unclassified.
- Coverage screening is relative to a sample baseline, without a control cohort, mappability correction or statistical CNV calls. GC correction may attenuate events correlated with GC content. Explicit neutral baseline contigs may be needed when much of the genome is altered.
- Synthetic tests establish execution and regression behavior, not biological sensitivity/specificity.
- No distribution license was selected automatically. License/maintainer/citation decisions and hosted CI remain publication checklist items.

## Release boundaries

Only the current diagram is distributed. Earlier diagram versions and prompt files are retained locally. Private species configs, sample sheets, standalone imports, biological results, Conda caches and generated release archives are excluded from source control. The local candidate branch/tag and release archive do not publish anything remotely.
