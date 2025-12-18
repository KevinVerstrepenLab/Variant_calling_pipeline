rule make_call_bam:
    input:
        bam = f"{BAM}/{{sample}}.dedup.bam",
        bai = f"{BAM}/{{sample}}.dedup.bam.bai"
    output:
        bam = f"{BAM}/{{sample}}.call.bam",
        bai = f"{BAM}/{{sample}}.call.bam.bai"
    conda: f"{ENVS}/bwa_samtools.yaml"
    threads: 16
    resources:
        mem_mb = 4000
    log:
        f"{LOGS}/callbam/{{sample}}.log"
    params:
        ppflag = ("-f 0x2" if config.get("call_bam", {}).get("require_proper_pair", True) else ""),
        excl   = config.get("call_bam", {}).get("exclude_flags", 3844),
        mapq   = config.get("call_bam", {}).get("mapq", 30)
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LOGS}/callbam
        # Keep: primary, mapped, non-dup, non-QCfail; require proper pairs; MAPQ≥threshold
        samtools view -@ {threads} -b {params.ppflag} -F {params.excl} -q {params.mapq} {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """
