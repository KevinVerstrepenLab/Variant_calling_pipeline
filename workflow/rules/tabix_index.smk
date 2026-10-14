rule tabix_index:
    input:
        f"{GVCF}/{{sample}}.g.vcf.gz"
    output:
        f"{GVCF}/{{sample}}.g.vcf.gz.tbi"
    conda: str(ENVS / "htslib.yaml")
    threads: 1
    resources:
        mem_mb=500
    log:
        f"{LOGS}/tabix/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')"
        tabix -f -p vcf '{input}' > '{log}' 2>&1
        """


rule index_legacy_gvcf:
    input:
        f"{GVCF}/{{stem}}.g.vcf"
    output:
        f"{GVCF}/{{stem}}.g.vcf.idx"
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=2000
    log:
        f"{LOGS}/legacy_index/{{stem}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')"
        gatk IndexFeatureFile -I '{input}' > '{log}' 2>&1
        test -s '{output}'
        """
