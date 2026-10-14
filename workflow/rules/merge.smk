# Merge filtered interval shards into genome-wide VCFs.
rule concat_snps:
    input:
        vcfs=expand(f"{SNPS_DIR}/filtered_snps.{{chunk}}.vcf.gz", chunk=CHUNKS),
        tbis=expand(f"{SNPS_DIR}/filtered_snps.{{chunk}}.vcf.gz.tbi", chunk=CHUNKS)
    output:
        vcf=f"{SNPS_DIR}/filtered_snps.merged.vcf.gz",
        tbi=f"{SNPS_DIR}/filtered_snps.merged.vcf.gz.tbi"
    conda: str(ENVS / "vcfops.yaml")
    threads: config.get("merge_threads", 2)
    resources:
        mem_mb=config.get("merge_mem_mb", 4000)
    log:
        f"{LOGS}/merge/snps_concat.log"
    params:
        snps_dir=lambda w, output: str(Path(output[0]).parent)

    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.snps_dir} "$(dirname '{log}')"
        bcftools concat --threads {threads} -Oz -o '{output.vcf}' {input.vcfs} > '{log}' 2>&1
        tabix -f -p vcf '{output.vcf}' >> '{log}' 2>&1
        """


rule concat_indels:
    input:
        vcfs=expand(f"{INDELS_DIR}/filtered_indels.{{chunk}}.vcf.gz", chunk=CHUNKS),
        tbis=expand(f"{INDELS_DIR}/filtered_indels.{{chunk}}.vcf.gz.tbi", chunk=CHUNKS)
    output:
        vcf=f"{INDELS_DIR}/filtered_indels.merged.vcf.gz",
        tbi=f"{INDELS_DIR}/filtered_indels.merged.vcf.gz.tbi"
    conda: str(ENVS / "vcfops.yaml")
    threads: config.get("merge_threads", 2)
    resources:
        mem_mb=config.get("merge_mem_mb", 4000)
    log:
        f"{LOGS}/merge/indels_concat.log"
    params:
        indels_dir=lambda w, output: str(Path(output[0]).parent)

    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.indels_dir} "$(dirname '{log}')"
        bcftools concat --threads {threads} -Oz -o '{output.vcf}' {input.vcfs} > '{log}' 2>&1
        tabix -f -p vcf '{output.vcf}' >> '{log}' 2>&1
        """
