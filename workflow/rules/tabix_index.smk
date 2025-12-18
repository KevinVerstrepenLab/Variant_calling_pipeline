rule tabix_index:
    input: f"{GVCF}/{{sample}}.g.vcf.gz"
    output: f"{GVCF}/{{sample}}.g.vcf.gz.tbi"
    conda: f"{ENVS}/htslib.yaml"
    log: f"{LOGS}/tabix/{{sample}}.log"
    threads: 1
    resources: mem_mb=500
    shell: "tabix -p vcf {input} 2> {log}"
