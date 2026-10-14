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
        mem_mb=config.get("markdup_mem_mb", 16000)
    log:
        f"{LOGS}/markdup/{{sample}}.log"
    params:
        bam_dir=lambda w, output: str(Path(output[0]).parent),
        java_opts=lambda w, resources: f"-Xms1g -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.bam_dir} "$(dirname '{log}')"
        gatk --java-options '{params.java_opts}' MarkDuplicatesSpark           --spark-master local[{threads}]           -I '{input.bam}' -O '{output.bam}' -M '{output.metrics}'           --conf 'spark.executor.cores={threads}'           > '{log}' 2>&1
        samtools index -@ {threads} '{output.bam}' '{output.bai}' 2>> '{log}'
        """
