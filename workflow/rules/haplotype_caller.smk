rule haplotype_caller_gvcf:
    input:
        bam = rules.make_call_bam.output.bam,
        bai = rules.make_call_bam.output.bai,
        ref = REF,
        fai = rules.ref_fai.output,
        dict = rules.ref_dict.output
    output:
        gvcf = f"{GVCF}/{{sample}}.g.vcf.gz"
    conda: f"{ENVS}/gatk.yaml"
    log: f"{LOGS}/haplotypecaller/{{sample}}.log"
    threads: config["hc_threads"]
    resources:
        mem_mb = config["hc_mem_mb"]
    params:
        ploidy = config.get("gatk", {}).get("sample_ploidy", 2),
        min_bq = config.get("gatk", {}).get("min_base_qual", 20),
        duscb  = "--dont-use-soft-clipped-bases" if config.get("gatk", {}).get("dont_use_soft_clipped", True) else ""
    log:
        f"logs/haplotypecaller/{{sample}}.log"
    shell:
        r"""
        mkdir -p {GVCF} {LOGS}/haplotypecaller
        gatk HaplotypeCaller \
          -R {input.ref} -I {input.bam} \
          -O {output.gvcf} -ERC GVCF \
          --sample-ploidy {params.ploidy} \
          --min-base-quality-score {params.min_bq} \
          {params.duscb} \
          2> {log}
        """
