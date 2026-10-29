rule haplotype_caller_gvcf:
    input:
        bam=rules.make_call_bam.output.bam,
        bai=rules.make_call_bam.output.bai,
        ref=REF,
        fai=rules.ref_fai.output,
        dict=rules.ref_dict.output
    output:
        gvcf=f"{GVCF}/{{sample}}.g.vcf.gz"
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in CALL_SAMPLES) if CALL_SAMPLES else r"(?!)"
    conda: str(ENVS / "gatk.yaml")
    threads: config.get("hc_threads", 1)
    resources:
        mem_mb=config.get("hc_mem_mb", 16000)
    log:
        f"{LOGS}/haplotypecaller/{{sample}}.log"
    params:
        gvcf_dir=lambda w, output: str(Path(output[0]).parent),
        staging_template=lambda w, output: str(Path(output[0]).parent / f".{w.sample}.hc.XXXXXX"),
        tmpdir=TMP_DIR,
        ploidy=config.get("gatk", {}).get("sample_ploidy", 2),
        min_bq=config.get("gatk", {}).get("min_base_qual", 20),
        duscb="--dont-use-soft-clipped-bases" if config.get("gatk", {}).get("dont_use_soft_clipped", True) else "",
        java_opts=lambda w, resources: f"-Xms1g -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.gvcf_dir:q} {params.tmpdir:q} "$(dirname {log:q})"
        staging_dir=$(mktemp -d {params.staging_template:q})
        trap 'rm -rf -- "$staging_dir"' EXIT
        staged_gvcf="$staging_dir/calls.g.vcf.gz"
        gatk --java-options {params.java_opts:q} HaplotypeCaller \
          -R {input.ref:q} -I {input.bam:q} \
          -O "$staged_gvcf" -ERC GVCF --native-pair-hmm-threads {threads} \
          --sample-ploidy {params.ploidy} \
          --min-base-quality-score {params.min_bq} \
          --tmp-dir {params.tmpdir:q} \
          {params.duscb} > {log:q} 2>&1
        # Only expose a reusable GVCF after GATK has finished writing both files.
        test -s "$staged_gvcf"
        test -s "$staged_gvcf.tbi"
        mv -- "$staged_gvcf.tbi" {output.gvcf:q}.tbi
        mv -- "$staged_gvcf" {output.gvcf:q}
        """
