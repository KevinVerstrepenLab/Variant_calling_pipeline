rule genotype_gvcfs:
    input:
        db   = rules.genomicsdbimport.output.db,
        ref  = REF,
        fai  = rules.ref_fai.output,
        dict = rules.ref_dict.output
    output:
        vcf = f"{JOINT}/Genotype.{{chunk}}.vcf.gz",
        tbi = f"{JOINT}/Genotype.{{chunk}}.vcf.gz.tbi"
    params:
        interval = lambda w: INTERVAL_BY_TAG[w.chunk],
        tmpdir = TMP_DIR,
        ploidy = config.get("gatk", {}).get("sample_ploidy", 2)
    conda: f"{ENVS}/gatk.yaml"
    log: f"{LOGS}/genotype/{{chunk}}.log"
    threads: config.get("hc_threads", 8)
    resources:
        mem_mb = config.get("hc_mem_mb", 16000)
    shell:
        r"""
        gatk GenotypeGVCFs \
          -R '{input.ref}' \
          -V 'gendb://{input.db}' \
          -L '{params.interval}' \
          --sample-ploidy {params.ploidy} \
          -O '{output.vcf}' \
          --tmp-dir {params.tmpdir} 2> {log}
        """
