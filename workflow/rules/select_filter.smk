# Select SNPs
rule select_snps:
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        ref = REF
    output:
        vcf = f"{SNPS_DIR}/{{chunk}}.snps.vcf.gz"
    conda: f"{ENVS}/gatk.yaml"
    threads: 1
    log: f"{LOGS}/select_snps/{{chunk}}.log"
    shell:
        r"""
        mkdir -p {SNPS_DIR} {LOGS}/select_snps
        gatk SelectVariants -V {input.vcf} -R {input.ref} -select-type SNP -O {output.vcf} 2> {log}
        """

# Filter SNPs (GUARDED JEXL)
rule filter_snps:
    input:
        vcf = rules.select_snps.output.vcf,
        ref = REF
    output:
        vcf = f"{SNPS_DIR}/filtered_snps.{{chunk}}.vcf.gz"
    params:
        # guarded version of: QD < 5 || FS > 30 || MQ < 50 || MQRankSum < -8 || ReadPosRankSum < -5 || SOR > 2
        expr = (
            '(vc.hasAttribute("QD") && QD < 5.0) || '
            '(vc.hasAttribute("FS") && FS > 30.0) || '
            '(vc.hasAttribute("MQ") && MQ < 50.0) || '
            '(vc.hasAttribute("MQRankSum") && MQRankSum < -8.0) || '
            '(vc.hasAttribute("ReadPosRankSum") && ReadPosRankSum < -5.0) || '
            '(vc.hasAttribute("SOR") && SOR > 2.0)'
        ),
        tmpdir = TMP_DIR
    conda: f"{ENVS}/gatk.yaml"
    threads: 1
    log: f"{LOGS}/filter_snps/{{chunk}}.log"
    shell:
        r"""
        mkdir -p {LOGS}/filter_snps
        gatk VariantFiltration \
          -V {input.vcf} -R {input.ref} \
          --filter-expression '{params.expr}' --filter-name snp_filter \
          -O {output.vcf} --tmp-dir {params.tmpdir} 2> {log}
        """

# Select INDELs
rule select_indels:
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        ref = REF
    output:
        vcf = f"{INDELS_DIR}/{{chunk}}.indels.vcf.gz"
    conda: f"{ENVS}/gatk.yaml"
    threads: 1
    log: f"{LOGS}/select_indels/{{chunk}}.log"
    shell:
        r"""
        mkdir -p {INDELS_DIR} {LOGS}/select_indels
        gatk SelectVariants -V {input.vcf} -R {input.ref} -select-type INDEL -O {output.vcf} 2> {log}
        """

# Filter INDELs (GUARDED JEXL)
rule filter_indels:
    input:
        vcf = rules.select_indels.output.vcf,
        ref = REF
    output:
        vcf = f"{INDELS_DIR}/filtered_indels.{{chunk}}.vcf.gz"
    params:
        # guarded version of: QD < 5 || FS > 200 || ReadPosRankSum < -10 || SOR > 10
        expr = (
            '(vc.hasAttribute("QD") && QD < 5.0) || '
            '(vc.hasAttribute("FS") && FS > 200.0) || '
            '(vc.hasAttribute("ReadPosRankSum") && ReadPosRankSum < -10.0) || '
            '(vc.hasAttribute("SOR") && SOR > 10.0)'
        ),
        tmpdir = TMP_DIR
    conda: f"{ENVS}/gatk.yaml"
    threads: 1
    log: f"{LOGS}/filter_indels/{{chunk}}.log"
    shell:
        r"""
        mkdir -p {LOGS}/filter_indels
        gatk VariantFiltration \
          -V {input.vcf} -R {input.ref} \
          --filter-expression '{params.expr}' --filter-name indel_filter \
          -O {output.vcf} --tmp-dir {params.tmpdir} 2> {log}
        """
