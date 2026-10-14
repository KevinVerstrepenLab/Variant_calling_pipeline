rule make_call_bam:
    input:
        bam=f"{BAM}/{{sample}}.dedup.bam",
        bai=f"{BAM}/{{sample}}.dedup.bam.bai"
    output:
        bam=f"{BAM}/{{sample}}.call.bam",
        bai=f"{BAM}/{{sample}}.call.bam.bai"
    conda: str(ENVS / "bwa_samtools.yaml")
    threads: config.get("call_bam_threads", config.get("bwa_threads", 6))
    resources:
        mem_mb=config.get("call_bam_mem_mb", 4000)
    log:
        f"{LOGS}/callbam/{{sample}}.log"
    params:
        bam_dir=lambda w, output: str(Path(output[0]).parent),
        ppflag=lambda w: "-f 0x2" if (
            LAYOUT[w.sample] == "paired" and
            config.get("call_bam", {}).get("require_proper_pair", True)
        ) else "",
        excl=config.get("call_bam", {}).get("exclude_flags", 3844),
        mapq=config.get("call_bam", {}).get("mapq", 30)
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.bam_dir} "$(dirname '{log}')"
        samtools view -@ {threads} -b {params.ppflag} -F {params.excl} -q {params.mapq} '{input.bam}' 2> '{log}'           | samtools sort -@ {threads} -o '{output.bam}' - 2>> '{log}'
        samtools index -@ {threads} '{output.bam}' '{output.bai}' 2>> '{log}'
        """
