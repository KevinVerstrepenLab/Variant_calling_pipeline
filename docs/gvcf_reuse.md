# Reusing existing GVCFs when adding samples

The workflow scans `<outdir>/gvcf` when the DAG is created. The public example
configuration enables `gvcf_reuse.include_existing_samples: true`, combining
every existing GVCF with samples in the input sheet. Set it to `false` to
restrict the cohort to sheet samples. Existing GVCFs bypass trimming,
alignment, duplicate marking and HaplotypeCaller.

## Safety condition: identical reference

GVCFs must never be combined when they were produced against different
reference sequences. The workflow therefore creates:

- `<outdir>/provenance/reference.sha256`
- `<outdir>/provenance/gvcf/<sample>.meta.json`
- `<outdir>/provenance/cohort.validated.json`

Before a GVCF enters `samples.map`, its contig names and lengths are checked
against the configured FASTA index and its reference identity is validated.
The preferred validation hashes the reference named in the GVCF `##reference`
header or HaplotypeCaller command-line header and compares it with the SHA-256 of the configured `ref`.

For older GVCFs whose header reference path is no longer accessible, strict
validation stops with an explanatory error. After independently confirming
that those legacy GVCFs used the identical reference, they can be adopted once
with:

```yaml
gvcf_reuse:
  allow_legacy_header_match: true
```

This fallback checks the complete ordered contig-name and contig-length list,
but cannot prove sequence identity. After metadata sidecars have been created,
set the option back to `false`.

## Adding samples

1. Keep the same `outdir` and exactly the same reference FASTA.
2. Retain all cohort samples in the sheet and add the new samples. If the sheet
   lists only new samples, set `gvcf_reuse.include_existing_samples: true` to
   retain the existing cohort.
3. Run Snakemake normally.

Existing indexed GVCFs are reused. Both compressed `.g.vcf.gz` (with `.tbi`)
and uncompressed `.g.vcf` (with `.idx`) are supported. Legacy
`<sample>.snps.indels.g.vcf` files are read directly without recompression.
Missing indexes are generated. New
samples proceed through FASTQ-to-GVCF processing. The sample map changes, so
GenomicsDB workspaces and all joint/downstream VCFs are rebuilt automatically.
Per-sample existing GVCFs are not regenerated.

## Important operational notes

- Filenames must match the actual sample header: `<sample>.g.vcf`,
  `<sample>.g.vcf.gz`, or their `<sample>.snps.indels.g.vcf[.gz]` equivalents.
- Multiple files identifying the same sample are rejected.
- Do not place unrelated GVCFs in `<outdir>/gvcf`.
- A zero-byte GVCF is rejected.
- If an existing GVCF is incomplete, remove it and its `.tbi` before rerunning.
- Adding samples necessarily reruns joint genotyping and downstream filtering;
  only the per-sample processing is reused.

Reused GVCFs are external DAG inputs, so upstream code/parameter changes do not
recall them. To recall samples with new calling settings, use a fresh output directory
or disable reuse. Validation checks the actual sample header and full gzip integrity.
Sidecars bind the reference fingerprint to a SHA-256 of the GVCF; a previously
validated sidecar can be used when the historical reference path disappears.

For large cohorts, `joint.batch_size: 50` limits readers open simultaneously.
`joint.consolidate: true` consolidates imported database batches. Provenance
validation runs as one cohort job with bounded parallel readers, using
`gvcf_reuse.validation_threads` (default 8). Snakemake tracks one atomic cohort
manifest, avoiding duplicated bookkeeping for every sidecar. The sample map
must match that manifest exactly before joint calling begins.
