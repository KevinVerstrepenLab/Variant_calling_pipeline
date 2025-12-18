# Merge filtered shards into genome-wide VCFs (bcftools)
rule concat_snps:
    input:
        expand(f"{SNPS_DIR}/filtered_snps.{{chunk}}.vcf.gz", chunk=CHUNKS)
    output:
        vcf = f"{SNPS_DIR}/filtered_snps.merged.vcf.gz",
        tbi = f"{SNPS_DIR}/filtered_snps.merged.vcf.gz.tbi"
    conda: f"{ENVS}/bwa_samtools.yaml"
    log: f"{LOGS}/merge/snps_concat.log"
    threads: 2
    shell:
        r"""
        bcftools concat -Oz -o {output.vcf} {input}
        tabix -p vcf {output.vcf} 2>> {log}
        """

rule concat_indels:
    input:
        expand(f"{INDELS_DIR}/filtered_indels.{{chunk}}.vcf.gz", chunk=CHUNKS)
    output:
        vcf = f"{INDELS_DIR}/filtered_indels.merged.vcf.gz",
        tbi = f"{INDELS_DIR}/filtered_indels.merged.vcf.gz.tbi"
    conda: f"{ENVS}/bwa_samtools.yaml"
    log: f"{LOGS}/merge/indels_concat.log"
    threads: 2
    shell:
        r"""
        bcftools concat -Oz -o {output.vcf} {input}
        tabix -p vcf {output.vcf} 2>> {log}
        """
