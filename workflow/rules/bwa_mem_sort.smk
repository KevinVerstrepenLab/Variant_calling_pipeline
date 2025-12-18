rule bwa_mem_sort:
    input:
        ref=REF,
        fai=rules.ref_fai.output,
        dict=rules.ref_dict.output,
        idx=rules.ref_bwa_index.output,
        r1=f"{TRIM}/{{sample}}_R1.paired.fq.gz",
        r2=f"{TRIM}/{{sample}}_R2.paired.fq.gz"
    output:
        bam=f"{BAM}/{{sample}}.sorted.bam"
    conda: f"{ENVS}/bwa_samtools.yaml"
    log: f"{LOGS}/markdup/{{sample}}.log"
    threads: config["bwa_threads"]
    resources: mem_mb=config["bwa_mem_mb"]
    log: f"logs/bwa/{{sample}}.log"
    shell: r"""
        mkdir -p {BAM} {LOGS}/bwa
        RG='@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA\tLB:lib1\tPU:unit1'
        bwa mem -t {threads} -R "$RG" {input.ref} {input.r1} {input.r2} 2> {log} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam} 2>> {log}
    """
