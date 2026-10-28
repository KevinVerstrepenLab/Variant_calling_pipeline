# Unified configuration and population analyses

`config/config.yaml` is the single run configuration. All relative paths resolve from the directory where Snakemake is run (normally the repository root), including metadata and sample-keep files. The workflow takes its analysis input directly from `<project.outdir>/final/filtered_snps.final.vcf.gz`; do not copy the standalone bundle's `datasets` or `output_root` sections into the run config.

## Configuration sections

| Section | Purpose |
| --- | --- |
| `project` | Output directory |
| `inputs` | Reference and sample sheet, or FASTQ discovery settings |
| `resources` | Threads and memory per processing stage |
| `read_processing` | fastp and calling-BAM filtering |
| `variant_calling` | GATK, joint genotyping and existing GVCF reuse |
| `filtering` | Variant hard filters, genotype masking and final missingness |
| `postprocessing` | Shared sample metadata, PCA/tree and LEA settings |
| `coverage_plots` | Optional BAM coverage screening and figure settings |

The bundled default leaves post-processing disabled. Enable `postprocessing.enabled`, then independently enable or disable `pca.enabled`, `phylogenetic_tree.enabled` and `lea.enabled`. For FASTQ discovery instead of a sheet, set `inputs.samples_tsv: null` and configure `inputs.fastq_dir`/`read_layout`. Existing flat configs continue to work. Legacy keys override nested equivalents when both are supplied; avoid mixing formats in a new run config.

Production code lives under `workflow/scripts/`; the integrated workflow does not read standalone-bundle YAML files. Standalone archives and species-specific local configs are not distributed in the release. Existing local configs can adopt the new sectioned format when needed.

## Migrating an existing configuration

| Previous key | Unified key |
| --- | --- |
| `ref` | `inputs.reference` |
| `samples_tsv`, `fastq_dir`, `read_layout`, `single_end_pattern` | Same names under `inputs` |
| `outdir` | `project.outdir` |
| `<stage>_threads`, `<stage>_mem_mb` | `resources.<stage>.threads`, `resources.<stage>.mem_mb` |
| `fastp_min_len`, `fastp_cut_front`, `fastp_cut_tail`, `fastp_cut_mean_quality` | Same names without `fastp_`, under `read_processing.fastp` |
| `call_bam` | `read_processing.call_bam` |
| `gatk`, `joint`, `gvcf_reuse` | Same names under `variant_calling` |
| `filters`, `ad_filter`, `final` | `filtering.hard`, `filtering.genotypes`, `filtering.final` |
| Bundle `pca_tree` | `postprocessing.variant_selection` for shared filters, `postprocessing.pca` for PCA, and `postprocessing.phylogenetic_tree` for tree options |
| Bundle `lea` | `postprocessing.lea` |
| Bundle rename/clade/palette paths | `postprocessing.metadata.rename`, `.clades`, `.clade_colors` |

The bundle's `filter` step is covered by existing workflow final filtering. PCA `geno` is the maximum missing fraction per site. LEA `max_missing` and `filtering.final.max_missing` are the minimum non-missing fraction per site. These thresholds have different directions.

## Sample metadata

Configure the following optional two-column TSV files once under `postprocessing.metadata`:

- `rename`: old sample ID and new label, without a header.
- `clades`: original ID or renamed label and clade name. A `sample<TAB>clade` header is accepted.
- `clade_colors`: clade name and colour such as `#336699`. A `clade<TAB>color` header is accepted.

Null values disable the metadata. Configured files must exist and be nonempty. For a branch-specific override, put the path under `postprocessing.variant_selection` (PCA/tree) or `postprocessing.lea`. Automatic colours use the same helper and full clade-name set in both stages. An explicit palette pins colours across runs with different clade sets. LEA ancestry-cluster colours are separate from clade annotation colours.

## PCA and trees

`postprocessing.variant_selection` controls shared MAF, missingness, sample selection and LD pruning. Separate `postprocessing.pca` and `postprocessing.phylogenetic_tree` sections control each analysis. The older combined `pca_tree` section remains supported.

PCA/tree preparation applies MAF, genotype missingness, SNP/biallelic and optional sample-keep filters, assigns unique variant IDs, and optionally prunes LD. PCA and tree generation can be enabled independently. For a PCA-only run, set `postprocessing.phylogenetic_tree.enabled: false`; PCA figures are still generated. PCA requires at least two components. PLINK2 LD pruning requires a sufficiently large sample set; disable `postprocessing.variant_selection.ld_prune.enabled` for small exploratory datasets if appropriate.

PCA component count and point labels are configured with `pca.components`, `pca.label_points` and `pca.label_size`. Tree tip labels use `phylogenetic_tree.label_tips` and `phylogenetic_tree.label_size` independently.

To display a subset, set `sample_list` independently in the PCA and TREE sections:

```yaml
postprocessing:
  pca:
    sample_list: metadata/pca_samples.txt
  phylogenetic_tree:
    sample_list: metadata/tree_samples.txt
```

Each file contains one sample ID per line, without a header. Original IDs and renamed labels are accepted; blank lines, duplicate IDs and lines starting with `#` are ignored. Paths resolve from the workflow working directory. Use `null` to display all samples. The lists may differ between PCA and TREE. Unknown IDs produce a warning; a list with no matching samples fails. TREE requires at least two matching tips.

These lists select the samples displayed after analysis. PCA coordinates and variance percentages are computed using the complete cohort selected by `variant_selection`, then the scores table and figures include only the requested PCA samples. The tree is inferred from that complete cohort and pruned for the figures. The `.nwk` file retains the full inferred tree; `.tree_info.tsv` records the displayed tips and the list used. To restrict the cohort before computing either analysis, use `variant_selection.keep` instead (PLINK keep-file format).

Editing either list schedules plotting and preserves existing PCA, distance and LEA results.

Tree methods are NJ and UPGMA on PLINK 1-IBS distances. The tree layout controls the PDF/PNG filename. `branch_length` controls the visualisation; the Newick file retains inferred lengths. `outgroup` refers to a tip name after renaming.

Outputs under `postprocessing/pca_tree/` include:

- `<prefix>.pca.eigenvec` and `.eigenval`;
- `<prefix>.tree.pca_scores.tsv`, `.tree.pca.pdf` and `.tree.pca.png`;
- `<prefix>.tree.nwk`, `.tree.tree_info.tsv` and `.tree-<layout>.pdf/png`;
- `<prefix>.tree.clade_colors.tsv` when clades are supplied.

## LEA

LEA preparation selects biallelic SNPs, applies its MAF and non-missingness thresholds, prunes LD, and exports PED/MAP and a pruned VCF. `run_snmf: false` stops after preparation.

sNMF tests `k_min` through `k_max` with the requested repetitions. `best_k: auto` selects the lowest mean cross-entropy; an integer selects that K from the tested range. The plotting range is intersected with the tested range. Empty intersections and invalid final K values fail during configuration.

Outputs under `postprocessing/lea/` include:

- `<prefix>.LEA.final.ped/.map` and `<prefix>.pruned.vcf.gz` with CSI index;
- `<prefix>.LEA.geno`, `.LEA.snmfProject` and the `.LEA.snmf/` model directory;
- `<prefix>.LEA.cross_entropy.tsv`, `.LEA.cross_entropy_summary.tsv` and `.LEA.cross_entropy.pdf`;
- `<prefix>.LEA.bestK.txt`, `.LEA.best.Q.tsv` and `.LEA.best.admixture.pdf`;
- K-specific Q/plot files and `<prefix>.LEA.K<min>_to_K<max>.admixture.pdf`;
- clade palettes and sample/clade assignments when metadata is supplied.

The stable `best` filenames let Snakemake track auto-selected K products directly.

## Scheduling, reuse and logs

Snakemake tracks the final SNP VCF/index, metadata files, scripts, parameters and main analysis products, including the fitted Q/G matrices. Five rules separate preparation and plotting for PCA/tree, and preparation, model fitting and plotting for LEA. Plot settings and sample annotations are excluded from preparation/model job parameters, so changing them schedules plotting alone. Each job uses `resources.postprocessing` thread/memory budgets; PLINK receives these limits explicitly. Each branch also retains internal checkpoints and intermediate files. Keep the `.LEA.snmf/` directory with the project.

Changing `best_k`, ordering or plot options reuses a matching sNMF model. Changing genotypes, K range, repetitions or seed invalidates the model. Models with missing or mismatching checkpoint provenance are recomputed. Shared palette helper changes invalidate plotting checkpoints.

To force computation inside a scheduled phase, set `postprocessing.force: true` and request its rule with `--forcerun`, for example `--forcerun postprocess_lea_model`. Reset `force` afterwards. A Snakemake force request alone still permits internal checkpoint reuse where intermediate outputs remain available. Plotting alone never refits a model.

Snakemake logs are under `logs/postprocessing/`. Each analysis directory also contains `commands.log` and `.checkpoints/`. Dependencies are provided by `workflow/envs/postprocessing.yaml` when running with `--use-conda`.
