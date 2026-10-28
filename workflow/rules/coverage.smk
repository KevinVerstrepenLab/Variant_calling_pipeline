import json
from configuration import coverage_settings

try:
    COVERAGE = coverage_settings(config)
except (ValueError, TypeError) as exc:
    raise WorkflowError(str(exc)) from exc
COVERAGE_DIR = f"{OUT}/CoveragePlots"
COVERAGE_INPUT_DIR = Path(COVERAGE["bam_dir"] or BAM)
COVERAGE_BAMS = {}
COVERAGE_EXISTING = set()
if COVERAGE["enabled"]:
    COVERAGE_BAMS = {path.name.removesuffix(".sorted.bam"): str(path.resolve())
                     for path in sorted(COVERAGE_INPUT_DIR.glob("*.sorted.bam"))}
    COVERAGE_EXISTING = set(COVERAGE_BAMS)
    if COVERAGE["bam_dir"] is None:
        for sample in CALL_SAMPLES:
            COVERAGE_BAMS.setdefault(sample, f"{BAM}/{sample}.sorted.bam")
COVERAGE_SAMPLES = sorted(COVERAGE_BAMS)
if COVERAGE["enabled"] and not COVERAGE_SAMPLES:
    raise WorkflowError("Coverage screening is enabled but no sorted BAMs or newly called samples are available")
COVERAGE_TARGETS = [f"{COVERAGE_DIR}/{sample}{suffix}" for sample in COVERAGE_SAMPLES
                    for suffix in (".png", ".pdf", ".windows.tsv.gz", ".candidate_regions.tsv", ".qc.json")]


rule coverage_plots:
    input:
        COVERAGE_TARGETS


if COVERAGE["enabled"] and COVERAGE_SAMPLES:
    rule coverage_reference_windows:
        input:
            ref=REF,
            fai=REF_FAI,
            script=str(WORKFLOW / "scripts" / "coverage_windows.py")
        output:
            bed=f"{COVERAGE_DIR}/reference/windows.bed",
            metadata=f"{COVERAGE_DIR}/reference/windows.tsv"
        params:
            settings=json.dumps({key: COVERAGE[key] for key in
                                 ("window_bp", "min_contig_bp", "exclude_contigs", "max_n_fraction")}, sort_keys=True)
        threads: 1
        resources:
            mem_mb=config.get("coverage_plots_mem_mb", 4000)
        conda: str(ENVS / "coverage.yaml")
        log: f"{LOGS}/coverage/reference_windows.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {log:q})"
            python {input.script:q} --reference {input.ref:q} --bed {output.bed:q} \
                --metadata {output.metadata:q} --settings {params.settings:q} > {log:q} 2>&1
            """

    rule coverage_depth:
        input:
            # Existing BAMs are external evidence: do not expand their alignment
            # producers or require historical FASTQs. Their stat signature below
            # triggers recomputation when an external BAM changes.
            bam=lambda w: [] if w.sample in COVERAGE_EXISTING else [COVERAGE_BAMS[w.sample]],
            bed=rules.coverage_reference_windows.output.bed
        output:
            regions=f"{COVERAGE_DIR}/depth/{{sample}}.regions.bed.gz",
            summary=f"{COVERAGE_DIR}/depth/{{sample}}.mosdepth.summary.txt",
            staged_bam=temp(f"{COVERAGE_DIR}/depth/{{sample}}.bam"),
            staged_index=temp(f"{COVERAGE_DIR}/depth/{{sample}}.bam.bai")
        params:
            prefix=lambda w, output: str(output.regions).removesuffix(".regions.bed.gz"),
            bam=lambda w: COVERAGE_BAMS[w.sample],
            source_signature=lambda w: (
                (Path(COVERAGE_BAMS[w.sample]).stat().st_size,
                 Path(COVERAGE_BAMS[w.sample]).stat().st_mtime_ns)
                if w.sample in COVERAGE_EXISTING else None
            ),
            mapq=COVERAGE["mapq"],
            flags=COVERAGE["exclude_flags"]
        threads: config.get("coverage_plots_threads", 4)
        resources:
            mem_mb=config.get("coverage_plots_mem_mb", 4000)
        conda: str(ENVS / "coverage.yaml")
        log: f"{LOGS}/coverage/{{sample}}.depth.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {output.regions:q})" "$(dirname {log:q})"
            ln -sfn {params.bam:q} {output.staged_bam:q}
            samtools index -@ {threads} -o {output.staged_index:q} {output.staged_bam:q} > {log:q} 2>&1
            mosdepth --threads {threads} --no-per-base --by {input.bed:q} \
                --mapq {params.mapq} --flag {params.flags} \
                {params.prefix:q} {output.staged_bam:q} >> {log:q} 2>&1
            """

    rule coverage_plot:
        input:
            regions=rules.coverage_depth.output.regions,
            metadata=rules.coverage_reference_windows.output.metadata,
            script=str(WORKFLOW / "scripts" / "plot_coverage.py")
        output:
            png=f"{COVERAGE_DIR}/{{sample}}.png",
            pdf=f"{COVERAGE_DIR}/{{sample}}.pdf",
            windows=f"{COVERAGE_DIR}/{{sample}}.windows.tsv.gz",
            candidates=f"{COVERAGE_DIR}/{{sample}}.candidate_regions.tsv",
            qc=f"{COVERAGE_DIR}/{{sample}}.qc.json"
        params:
            prefix=lambda w, output: str(Path(output.png).with_suffix("")),
            settings=json.dumps({key: value for key, value in COVERAGE.items()
                                 if key not in ("enabled", "bam_dir", "mapq", "exclude_flags")}, sort_keys=True)
        threads: 1
        resources:
            mem_mb=config.get("coverage_plots_mem_mb", 4000)
        conda: str(ENVS / "coverage.yaml")
        log: f"{LOGS}/coverage/{{sample}}.plot.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {log:q})"
            python {input.script:q} --regions {input.regions:q} --metadata {input.metadata:q} \
                --prefix {params.prefix:q} --sample {wildcards.sample:q} \
                --settings {params.settings:q} > {log:q} 2>&1
            """
