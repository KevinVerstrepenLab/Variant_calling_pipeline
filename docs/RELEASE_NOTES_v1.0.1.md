# v1.0.1 — Reliable output publication and Java job concurrency

Released 16 November 2026.

This patch improves interrupted-run recovery and adds a direct concurrency limit for GATK jobs. It includes the HaplotypeCaller and GVCF-reuse fixes made after the original v1.0.0 tag.

## Changes

- Every GATK rule requests one `java_jobs` slot. Set `--resources java_jobs=8` to allow at most eight Java jobs concurrently, subject to CPU and memory limits.
- The README includes complete commands for a 32-core workflow, explains why per-job memory can limit concurrency, and shows optional memory/import-thread overrides.
- HaplotypeCaller publishes reusable GVCFs only after nonempty calls and indexes are ready.
- Duplicate marking stages and validates the BAM, metrics and indexes before publishing the BAM. Failed attempts clean their staging directories, and existing Spark partitions do not block a retry.
- Reused GVCF metadata is ordered consistently, avoiding spurious cohort rebuilds.
- Coverage annotations and legends have improved placement.
- Release publication runs the regression suite, lint and dry-run checks before creating the GitHub release, source archives and SHA-256 checksums.

## Run or resume

From the repository root, after configuring your input paths and activating Snakemake:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml --use-conda \
  --cores 32 --resources mem_mb=128000 java_jobs=8 \
  --printshellcmds --latency-wait 30 --rerun-incomplete
```

Use this example only when 32 cores and 128,000 MB are available to the workflow after reserving capacity for the system and other workloads. Increase or decrease `java_jobs` with the available RAM and per-rule memory requests. It caps GATK processes, not JVM thread counts. See the [README](https://github.com/KevinVerstrepenLab/Variant_calling_pipeline/blob/v1.0.1/README.md#control-the-number-of-concurrent-java-jobs) for the tuning examples.

## Validation

All 55 regression cases passed across the installed test environments. Structural checks, Snakemake 9.5.1 lint and a clean-source mixed paired/single-end dry-run passed. A fresh synthetic paired/single-end run exercised trimming through final VCFs, including actual GATK duplicate marking and HaplotypeCaller. A large biological cohort using existing GVCFs completed calling/filtering and PCA/tree analysis; its configured dry-run reported no work remaining. Synthetic tests also exercised coverage and LEA analysis.

The Java slot limit is a scheduling control. Smaller per-job memory settings are optional and require sufficient memory for the actual cohort. These checks establish execution and regression behavior; they do not measure biological sensitivity/specificity or guarantee a speedup. Conda environments were tested as installed, without a fresh solve of every environment.

Private configurations, sample sheets, biological results, logs, caches and archived figures are excluded from the source assets. No new distribution license is assigned by this patch.
