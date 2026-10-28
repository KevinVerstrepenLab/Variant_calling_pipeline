# Local coverage screening from sorted BAMs

Enable `coverage_plots.enabled: true` in the same run configuration. The workflow writes each sample's figures and screening tables under `<outdir>/CoveragePlots/`. This branch starts from `.sorted.bam` files and does not require variants to be called first.

```bash
snakemake coverage_plots --snakefile workflow/Snakefile \
  --configfile config/config.yaml --use-conda --cores 8
```

With the default `bam_dir: null`, the analysis includes every existing `<outdir>/bam/*.sorted.bam` plus the sorted BAMs that will be produced for newly called samples. Existing BAMs include samples absent from the current sheet and samples with reused GVCFs. Existing GVCF samples without a sorted BAM are skipped; coverage screening does not realign them. An explicit `bam_dir` selects only existing `.sorted.bam` files in that directory. Discovery is not recursive. Re-run Snakemake after adding files to include them.

Existing BAMs are treated as external evidence: coverage jobs track their resolved path, size and nanosecond modification time as parameters, avoiding expansion of alignment jobs and requirements for old FASTQs. Changing an existing BAM's size or modification time schedules depth recomputation. Preserve accurate timestamps when replacing files. New BAMs retain their alignment dependencies. Temporary BAM symlinks and indexes are created in `CoveragePlots/depth/`, leaving existing source BAM directories unchanged.

## Method

[Mosdepth](https://github.com/brentp/mosdepth) calculates mean depth in fixed-size genomic windows directly from sorted BAMs. The workflow uses its full CIGAR/overlapping-mate handling, disables per-base output, and applies configurable mapping-quality and alignment-flag filters. Zero-depth windows are retained. Duplicate alignments are excluded only if they are already flagged: the requested `.sorted.bam` files precede duplicate marking in this pipeline, so unmarked PCR duplicates can still inflate local depth.

Reference windows record GC content and the fraction of ambiguous bases. Short/excluded contigs are omitted; ambiguous windows remain visible as raw coverage but are masked from normalization and candidate detection. The baseline is a window-length-weighted median depth, using all eligible contigs or the exact names in `baseline_contigs`. This avoids normalizing every chromosome independently and preserves chromosome-scale changes relative to the rest of the genome. It also avoids giving short chromosomes the same weight as long chromosomes.

Optional GC correction sorts positive-depth baseline windows by GC, groups them into up to `gc_bins` equal-count groups (at least ten windows per group), calculates median depth per group and interpolates a GC-depth trend. Each depth is divided by the trend relative to its median. Correction is skipped, with the reason recorded, when there are fewer than `min_gc_windows` training windows or less than 0.02 GC-fraction variation. This simple robust screen is not an implementation of CNVkit's full correction model. The QC file records whether correction was applied.

Normalized values are `log2((corrected_depth + pseudocount) / (corrected_baseline + pseudocount))`. No integer ploidy states are assigned, and high coverage values are not discarded. Raw depths and raw log2 ratios are retained alongside corrected values. The pseudocount makes zero coverage finite and slightly attenuates ratios at low depth.

A centered rolling median smooths the signal within each contig and each contiguous run of valid windows. It uses shorter windows at run ends and never crosses masked windows or contig boundaries. Consecutive smoothed windows above/below the configured threshold are reported as candidate regions if their count meets `min_region_windows`. Coordinates are 0-based, end-exclusive. These regions are a visual screening aid, without statistical significance or exact breakpoint estimates.

## Configuration and window size

Full defaults are in section 7 of `config/config.yaml`. The public template leaves coverage disabled; enable it in your local run config when needed.

```yaml
coverage_plots:
  enabled: true
  window_bp: 10000
  smooth_windows: 5
  min_region_windows: 3
  log2_threshold: 0.40
  gc_correct: true
  exclude_contigs: []
  baseline_contigs: []
  plot_min_contig_bp: null
  chromosomes_per_row: 16
  outlier_quantiles: [0.01, 0.99]
```

The 10 kb default is a conservative starting choice for exploratory yeast coverage plots, rather than a universal optimum. Try 5 kb for finer local detail with adequate depth, or 20 kb if the track is noisy. A 1 kb exploratory track can help resolve smaller events at high depth, but needs closer assessment of noise, mapping and repeats. As a reference point, [CNVkit's WGS documentation](https://cnvkit.readthedocs.io/en/stable/nonhybrid.html) discusses bins of at least 1 kb at 30x and larger bins for lower coverage; its recommendations accompany a more complete bias-correction/segmentation workflow.

At the defaults, smoothing spans approximately 50 kb and candidate regions require at least three 10 kb windows (30 kb). Smoothing can attenuate small events and move apparent boundaries. For finer resolution, adjust `window_bp`, `smooth_windows` and `min_region_windows` together; use `smooth_windows: 1` to inspect unsmoothed changes. Changing these plot/screen parameters reuses existing depth windows; changing `window_bp`, mapping filters or the selected reference regions recalculates depth.

List mitochondrial contigs and other inappropriate baseline regions by exact name in `exclude_contigs`. Use known neutral `baseline_contigs` when many chromosomes may be altered. A genome-wide baseline assumes most selected sequence is neutral. Uniform genome doubling cannot be identified by this normalization. Low baseline depth is flagged in the title and QC file; a zero median baseline stops with an explicit error.

Sequencing depth also reflects GC content, repeats, mappability and library artifacts. [CNVkit's bias documentation](https://cnvkit.readthedocs.io/en/stable/bias.html) describes these confounders. GC and genomic location can be correlated, so single-sample GC correction may attenuate true events: compare raw and corrected tracks or disable it. For formal CNV inference, a tool such as CNVkit with an appropriate control/reference, accessibility masks and segmentation provides a stronger follow-up. Candidate coverage differences here are not confirmed aneuploidy or global ploidy estimates.

## Readable figures: outliers, short sequences and wrapping

These settings affect figures only. They do not discard analysis windows, change the coverage baseline, alter candidate regions, or rerun mosdepth:

- `plot_min_contig_bp`: minimum sequence length displayed on the x-axis. `null` or an omitted value means **50,000 bp**; exactly 50 kb is included. Use `0` to display every analyzed sequence. This is independent of `min_contig_bp`, which excludes short sequences before coverage analysis. To include sequences shorter than that analysis cutoff, lower `min_contig_bp` as well.
- `chromosomes_per_row`: maximum displayed sequences per row, **16** by default. Sequences retain reference order after the display-size filter. Seventeen sequences use two rows; 33 use three. Each row has a raw-depth track and a log2 track with its own sequence x-axis. `height` is the height per sequence row, so the total image height increases with wrapping. Y limits are shared across rows for comparison, and cumulative Mb coordinates restart within each row.
- `outlier_quantiles`: robust y-range limits, **`[0.01, 0.99]`** by default. Central-range limits are calculated from displayed windows, using inward-rounded percentile ranks to exclude isolated extreme points. The raw track always includes the baseline, and the log2 track always includes zero and the screening thresholds. Out-of-range observations are replaced in the figure by upward/downward edge triangles and counted in annotations. Smoothed traces break outside the visible range. No observations are removed from the saved tables. Use `null` for the complete y-range, or adjust the two quantiles for more/less trimming. Broad changes can still set the range; these limits are a display aid, not an outlier classification or biological filter.

The QC JSON records the displayed and hidden contigs, number of sequence rows, y limits, and out-of-range counts for each track. If no sequence passes the display-length cutoff, plotting stops with a message to lower `plot_min_contig_bp`.

## Outputs and reuse

For each sample:

- `<sample>.png` and `<sample>.pdf`: raw mean depth above, raw/normalized log2 ratios and a smoothed trace below; displayed contigs are shown in reference order, wrapped at `chromosomes_per_row`.
- `<sample>.windows.tsv.gz`: window coordinates, reference GC/ambiguity metadata, raw/corrected depth, raw/corrected log2 ratios and smoothed ratios.
- `<sample>.candidate_regions.tsv`: candidate gain/loss regions; header only if none are detected.
- `<sample>.qc.json`: baseline depth, low-depth flag, correction status, window counts and parameters.

Shared reference metadata is in `CoveragePlots/reference/`; compressed depth windows and mosdepth summaries are in `CoveragePlots/depth/`. Logs are under `logs/coverage/`. Depth calculation and plotting are separate Snakemake jobs, so changing thresholds, smoothing or display options does not rerun mosdepth. Resources are configured with `resources.coverage_plots` (or legacy `coverage_plots_threads`/`coverage_plots_mem_mb`). Conda dependencies are in `workflow/envs/coverage.yaml`.
