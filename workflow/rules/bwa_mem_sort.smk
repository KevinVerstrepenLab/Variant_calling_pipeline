from common import trimmed_reads

rule bwa_mem_sort:
    input:
        ref=REF,
        fai=rules.ref_fai.output,
        dict=rules.ref_dict.output,
        idx=rules.ref_bwa_index.output,
        reads=lambda wc: trimmed_reads(wc.sample, LAYOUT, TRIM)
    output:
        bam=f"{BAM}/{{sample}}.sorted.bam",
        bai=f"{BAM}/{{sample}}.sorted.bam.bai"
    conda:
        str(ENVS / "bwa_samtools.yaml")
    threads:
        config.get("bwa_threads", 8)
    resources:
        mem_mb=config.get("bwa_mem_mb", 12000)
    log:
        f"{LOGS}/bwa/{{sample}}.log"
    params:
        bam_dir=lambda w, output: str(Path(output[0]).parent),
        reads=lambda wc, input: " ".join(
            f"'{read}'" for read in input.reads
        ),
        rg=lambda wc: (
            f"@RG\\tID:{wc.sample}"
            f"\\tSM:{wc.sample}"
            f"\\tPL:ILLUMINA"
            f"\\tLB:lib1"
            f"\\tPU:unit1"
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.bam_dir}" "$(dirname "{log}")"

        bwa mem \
          -t {threads} \
          -R '{params.rg}' \
          '{input.ref}' \
          {params.reads} \
          2> '{log}' \
        | samtools sort \
            -@ {threads} \
            -o '{output.bam}' \
            - \
            2>> '{log}'

        samtools quickcheck -v '{output.bam}' 2>> '{log}'

        samtools index \
          -@ {threads} \
          '{output.bam}' \
          '{output.bai}' \
          2>> '{log}'
        """
