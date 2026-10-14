rule select_snps:
    input:
        vcf=rules.genotype_gvcfs.output.vcf,
        tbi=rules.genotype_gvcfs.output.tbi,
        ref=REF
    output:
        vcf=f"{SNPS_DIR}/{{chunk}}.snps.vcf.gz",
        tbi=f"{SNPS_DIR}/{{chunk}}.snps.vcf.gz.tbi"
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=config.get("select_mem_mb", 4000)
    log:
        f"{LOGS}/select_snps/{{chunk}}.log"
    params:
        snps_dir=lambda w, output: str(Path(output[0]).parent),
        java_opts=lambda w, resources: f"-Xms512m -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.snps_dir} "$(dirname '{log}')"
        gatk --java-options '{params.java_opts}' SelectVariants           -V '{input.vcf}' -R '{input.ref}' -select-type SNP -O '{output.vcf}'           > '{log}' 2>&1
        test -s '{output.tbi}'
        """


rule filter_snps:
    input:
        vcf=rules.select_snps.output.vcf,
        tbi=rules.select_snps.output.tbi,
        ref=REF
    output:
        vcf=f"{SNPS_DIR}/filtered_snps.{{chunk}}.vcf.gz",
        tbi=f"{SNPS_DIR}/filtered_snps.{{chunk}}.vcf.gz.tbi"
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=config.get("filter_mem_mb", 4000)
    log:
        f"{LOGS}/filter_snps/{{chunk}}.log"
    params:
        expr=lambda w: config.get("filters", {}).get("snp", 'QD < 5.0 || FS > 30.0 || MQ < 50.0 || MQRankSum < -8.0 || ReadPosRankSum < -5.0 || SOR > 2.0'),
        tmpdir=TMP_DIR,
        java_opts=lambda w, resources: f"-Xms512m -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')" '{params.tmpdir}'
        gatk --java-options '{params.java_opts}' VariantFiltration           -V '{input.vcf}' -R '{input.ref}'           --filter-expression '{params.expr}' --filter-name snp_filter           -O '{output.vcf}' --tmp-dir '{params.tmpdir}'           > '{log}' 2>&1
        test -s '{output.tbi}'
        """


rule select_indels:
    input:
        vcf=rules.genotype_gvcfs.output.vcf,
        tbi=rules.genotype_gvcfs.output.tbi,
        ref=REF
    output:
        vcf=f"{INDELS_DIR}/{{chunk}}.indels.vcf.gz",
        tbi=f"{INDELS_DIR}/{{chunk}}.indels.vcf.gz.tbi"
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=config.get("select_mem_mb", 4000)
    log:
        f"{LOGS}/select_indels/{{chunk}}.log"
    params:
        indels_dir=lambda w, output: str(Path(output[0]).parent),
        java_opts=lambda w, resources: f"-Xms512m -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.indels_dir} "$(dirname '{log}')"
        gatk --java-options '{params.java_opts}' SelectVariants           -V '{input.vcf}' -R '{input.ref}' -select-type INDEL -O '{output.vcf}'           > '{log}' 2>&1
        test -s '{output.tbi}'
        """


rule filter_indels:
    input:
        vcf=rules.select_indels.output.vcf,
        tbi=rules.select_indels.output.tbi,
        ref=REF
    output:
        vcf=f"{INDELS_DIR}/filtered_indels.{{chunk}}.vcf.gz",
        tbi=f"{INDELS_DIR}/filtered_indels.{{chunk}}.vcf.gz.tbi"
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=config.get("filter_mem_mb", 4000)
    log:
        f"{LOGS}/filter_indels/{{chunk}}.log"
    params:
        expr=lambda w: config.get("filters", {}).get("indel", 'QD < 5.0 || FS > 200.0 || ReadPosRankSum < -10.0 || SOR > 10.0'),
        tmpdir=TMP_DIR,
        java_opts=lambda w, resources: f"-Xms512m -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')" '{params.tmpdir}'
        gatk --java-options '{params.java_opts}' VariantFiltration           -V '{input.vcf}' -R '{input.ref}'           --filter-expression '{params.expr}' --filter-name indel_filter           -O '{output.vcf}' --tmp-dir '{params.tmpdir}'           > '{log}' 2>&1
        test -s '{output.tbi}'
        """
