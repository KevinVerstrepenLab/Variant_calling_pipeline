rule genotype_gvcfs:
    input:
        db=rules.genomicsdbimport.output.db,
        ref=REF,
        fai=rules.ref_fai.output,
        dict=rules.ref_dict.output
    output:
        vcf=f"{JOINT}/Genotype.{{chunk}}.vcf.gz",
        tbi=f"{JOINT}/Genotype.{{chunk}}.vcf.gz.tbi"
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=config.get("genotype_mem_mb", config.get("hc_mem_mb", 16000))
    log:
        f"{LOGS}/genotype/{{chunk}}.log"
    params:
        joint_dir=lambda w, output: str(Path(output[0]).parent),
        interval=lambda w: INTERVAL_BY_TAG[w.chunk],
        tmpdir=TMP_DIR,
        ploidy=config.get("gatk", {}).get("sample_ploidy", 2),
        java_opts=lambda w, resources: f"-Xms1g -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.joint_dir} "$(dirname '{log}')" '{params.tmpdir}'
        gatk --java-options '{params.java_opts}' GenotypeGVCFs           -R '{input.ref}'           -V 'gendb://{input.db}'           -L '{params.interval}'           --sample-ploidy {params.ploidy}           -O '{output.vcf}'           --tmp-dir '{params.tmpdir}'           > '{log}' 2>&1
        test -s '{output.tbi}'
        """
