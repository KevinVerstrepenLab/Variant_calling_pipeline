# Variant Calling Pipeline (Snakemake)

Call and jointly genotype variants from new paired-end or single-end reads together with existing per-sample GVCFs. Existing GVCFs are reused directly; only new samples need trimming, alignment and HaplotypeCaller.

![Variant calling workflow: new reads and existing GVCFs join a validated cohort before joint genotyping](VariantCalling2026.png)

## Configure a run

Run commands from the repository root. Copy `config/config.yaml` to your own configuration and set:

```yaml
ref: data/reference.fasta
samples_tsv: config/samples_examples.tsv
outdir: results
```

Replace the example paths with your actual files. Relative paths resolve from the working directory. The reference directory must be writable for FASTA, sequence-dictionary and BWA indexes.

The tab-separated sample sheet describes the new samples:

```text
sample      layout  R1                                  R2
isolate_01  paired  data/reads/isolate_01_R1.fastq.gz    data/reads/isolate_01_R2.fastq.gz
isolate_02  single  data/reads/isolate_02.fastq.gz
```

Use `config/samples_examples.tsv` as a template. `layout` accepts `paired` or `single`; leave `R2` empty for single-end reads. Older sheets containing `sample`, `R1`, and `R2` remain supported: layout is inferred from whether `R2` is populated. Sample names must be unique.

## Reuse an existing cohort

Keep existing GVCFs in `<outdir>/gvcf/`. Files can be regular files or symlinks. Supported names are:

- `<sample>.g.vcf.gz` with `<sample>.g.vcf.gz.tbi`;
- `<sample>.g.vcf` with `<sample>.g.vcf.idx`;
- `<sample>.snps.indels.g.vcf[.gz]` with the corresponding `.idx` or `.tbi`.

Uncompressed legacy GVCFs are read directly without recompression. Missing indexes are generated. Filenames must match the sample names in the GVCF headers; duplicate sample identities are rejected.

The public configuration enables:

```yaml
gvcf_reuse:
  enabled: true
  include_existing_samples: true
  validation_threads: 8
  allow_legacy_header_match: false
```

With these settings, the cohort includes every existing GVCF in `outdir/gvcf`, plus all samples in the input sheet. Existing samples do not need FASTQs in the sheet. If a sheet sample already has a GVCF, that GVCF is reused too. Set `include_existing_samples: false` to restrict the cohort to sheet samples. Use a fresh output directory or disable reuse to recall samples with changed calling settings.

Every selected GVCF is validated for sample identity, ordered contig names and lengths, and reference provenance. The workflow hashes the configured reference and records GVCF checksums. Compressed inputs are also checked for gzip integrity. Reference paths can come from either `##reference` or the GATK HaplotypeCaller command-line header.

If a historical reference path is unavailable, strict validation stops. After independently confirming that the legacy GVCFs used the identical reference sequence, enable `allow_legacy_header_match: true` for adoption. Matching contig names and lengths alone cannot prove sequence identity. See [GVCF reuse and provenance](docs/gvcf_reuse.md).

## Run

Install Snakemake and Conda, then inspect the planned jobs:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml --dry-run --cores 8
```

Execute the workflow:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda --cores 8 --resources mem_mb=32000 \
  --printshellcmds --latency-wait 30 --rerun-incomplete
```

Conda environments are defined under `workflow/envs`. Tune threads and memory in the configuration for your machine. Large cohorts use interval-based joint calling and bounded imports (`joint.batch_size`, default 50); `joint.consolidate` merges the database batches.

## Processing stages

1. **fastp:** trim adapters and low-quality bases; write read QC reports.
2. **BWA-MEM and samtools:** align, sort and index BAMs.
3. **GATK MarkDuplicatesSpark:** mark duplicate alignments.
4. **Calling BAM:** apply mapping-quality and flag filters. Proper-pair filtering applies only to paired-end libraries.
5. **GATK HaplotypeCaller:** generate GVCFs for new samples.
6. **Cohort validation:** combine existing and new GVCFs, validate provenance and build the complete sample map.
7. **GenomicsDBImport and GenotypeGVCFs:** import and jointly genotype each configured interval.
8. **Variant filtering:** select SNPs and INDELs, apply configurable hard filters, and concatenate interval VCFs.
9. **Final filtering:** optionally mask genotypes by allele depth, then apply sample and site missingness thresholds.

Calling-BAM and variant-filter defaults should be evaluated for the species, ploidy and divergence of your samples. See [filter validation guidance](docs/filter_validation.md). Allele-depth masking and heterozygote removal are disabled in the example configuration. AD masking assumes biallelic records; enable `ad_filter.split_multiallelic` when splitting is needed. `final.sample_max_missing` is a maximum missing fraction per sample, while `final.max_missing` is the minimum fraction of non-missing genotypes required at a site (vcftools semantics).

Adding samples rebuilds GenomicsDB and downstream joint results for the complete cohort. Reused per-sample GVCFs remain unchanged.

## Outputs

Generated results are under `outdir`:

- `trimmed/` and `bam/` for newly processed samples;
- `gvcf/` for existing and newly generated per-sample GVCFs and indexes;
- `provenance/reference.sha256`, per-sample sidecars in `provenance/gvcf/`, and `provenance/cohort.validated.json`;
- `joint/samples.map`, `genomicsdb/`, `joint/`, `snps/`, and `indels/`;
- `final/filtered_snps.final.vcf.gz` and `final/filtered_indels.final.vcf.gz`, each with a `.tbi` index;
- `logs/`, including cohort validation, import, genotyping and filtering logs.

SNP and INDEL missingness filters run independently, so the final sample sets can differ. The release includes portable example configuration and sample data paths; local datasets, private configurations, tests and runtime caches are excluded.
