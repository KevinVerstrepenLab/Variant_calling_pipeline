rule ref_fai:
    input:
        REF
    output:
        REF_FAI
    conda: str(ENVS / "bwa_samtools.yaml")
    threads: 1
    resources:
        mem_mb=500
    log:
        f"{LOGS}/ref/faidx.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')"
        samtools faidx '{input}' > '{log}' 2>&1
        """


rule ref_dict:
    input:
        REF
    output:
        REF_DICT
    conda: str(ENVS / "gatk.yaml")
    threads: 1
    resources:
        mem_mb=1000
    log:
        f"{LOGS}/ref/dict.log"
    params:
        java_opts=lambda w, resources: f"-Xms256m -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')"
        dict_tmpdir=$(mktemp -d "$(dirname '{output}')/.dict.XXXXXX")
        trap 'rm -rf "$dict_tmpdir"' EXIT
        gatk --java-options '{params.java_opts}' CreateSequenceDictionary \
          -R '{input}' -O "$dict_tmpdir/reference.dict" > '{log}' 2>&1
        mv "$dict_tmpdir/reference.dict" '{output}'
        """


rule ref_bwa_index:
    input:
        REF
    output:
        f"{REF}.amb",
        f"{REF}.ann",
        f"{REF}.bwt",
        f"{REF}.pac",
        f"{REF}.sa"
    conda: str(ENVS / "bwa_samtools.yaml")
    threads: 1
    resources:
        mem_mb=2000
    log:
        f"{LOGS}/ref/bwa_index.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')"
        bwa index '{input}' > '{log}' 2>&1
        """
