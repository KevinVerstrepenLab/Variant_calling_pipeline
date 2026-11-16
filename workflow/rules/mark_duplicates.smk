rule mark_duplicates:
    input:
        bam=f"{BAM}/{{sample}}.sorted.bam",
        bai=f"{BAM}/{{sample}}.sorted.bam.bai"
    output:
        bam=f"{BAM}/{{sample}}.dedup.bam",
        bai=f"{BAM}/{{sample}}.dedup.bam.bai",
        metrics=f"{BAM}/{{sample}}.dedup.metrics.txt"
    conda: str(ENVS / "gatk.yaml")
    threads: config.get("markdup_threads", 6)
    resources:
        java_jobs=1,
        mem_mb=config.get("markdup_mem_mb", 16000)
    log:
        f"{LOGS}/markdup/{{sample}}.log"
    params:
        bam_dir=lambda w, output: str(Path(output[0]).parent),
        staging_template=lambda w, output: str(Path(output[0]).parent / f".{w.sample}.markdup.XXXXXX"),
        tmpdir=TMP_DIR,
        java_opts=lambda w, resources: f"-Xms1g -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.bam_dir:q} {params.tmpdir:q} "$(dirname {log:q})"
        staging_dir=$(mktemp -d {params.staging_template:q})
        trap 'rm -rf -- "$staging_dir"' EXIT
        staged_bam="$staging_dir/marked.bam"
        staged_metrics="$staging_dir/metrics.txt"
        gatk --java-options {params.java_opts:q} MarkDuplicatesSpark \
          --spark-master local[{threads}] \
          -I {input.bam:q} -O "$staged_bam" -M "$staged_metrics" \
          --tmp-dir {params.tmpdir:q} \
          --conf 'spark.executor.cores={threads}' > {log:q} 2>&1
        samtools quickcheck -v "$staged_bam" 2>> {log:q}
        samtools index -@ {threads} "$staged_bam" "$staged_bam.bai" 2>> {log:q}
        # Publish the BAM last, after marking, validation and indexing succeed.
        test -s "$staged_bam"
        test -s "$staged_bam.bai"
        test -s "$staged_metrics"
        mv -- "$staged_metrics" {output.metrics:q}
        mv -- "$staged_bam.bai" {output.bai:q}
        if [ -f "$staged_bam.sbi" ]; then
            mv -- "$staged_bam.sbi" {output.bam:q}.sbi
        fi
        mv -- "$staged_bam" {output.bam:q}
        """
