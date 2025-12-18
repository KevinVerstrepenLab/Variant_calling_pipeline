rule mark_duplicates:
    input:
        f"{BAM}/{{sample}}.sorted.bam"
    output:
        bam = f"{BAM}/{{sample}}.dedup.bam",
        bai = f"{BAM}/{{sample}}.dedup.bam.bai",
        metrics = f"{BAM}/{{sample}}.dedup.metrics.txt"
    conda: f"{ENVS}/gatk.yaml"
    log: f"{LOGS}/markdup/{{sample}}.log"
    threads: config["markdup_threads"]
    resources:
        mem_mb = config["markdup_mem_mb"]
    log:
        f"logs/markdup/{{sample}}.log"
    shell:
        r"""
        mkdir -p logs/markdup
        gatk MarkDuplicatesSpark \
          --spark-master local[{threads}] \
          -I {input} -O {output.bam} -M {output.metrics} \
          --conf 'spark.executor.cores={threads}' \
          2> {log}
        samtools index {output.bam} 2>> {log}
        """
