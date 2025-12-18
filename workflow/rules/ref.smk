rule ref_fai:
    input: REF
    output: REF + ".fai"
    conda: f"{ENVS}/bwa_samtools.yaml"
    log:   f"{LOGS}/ref/faidx.log"
    threads: 1
    resources: mem_mb=500
    shell: "mkdir -p {LOGS}/ref; samtools faidx {input} 2> {log}"

rule ref_dict:
    input:
        REF
    output:
        REF_DICT
    conda:
        f"{ENVS}/gatk.yaml"
    log:
        f"{LOGS}/ref/dict.log"
    threads: 1
    resources:
        mem_mb=1000
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LOGS}/ref

        # Try with -O (works on most GATK 4 / Picard wrappers)
        if gatk CreateSequenceDictionary -R {input} -O {output} 2> {log}; then
            :
        else
            echo "[WARN] Retrying CreateSequenceDictionary without -O" >> {log}
            rm -f {output} || true
            gatk CreateSequenceDictionary -R {input} 2>> {log}
            # If it wrote next to the FASTA, link it to the expected output path
            if [ -e "{input}.dict" ]; then
              ln -sf "{input}.dict" "{output}"
            fi
        fi

        # Verify we have the expected file
        if [ ! -e "{output}" ]; then
          echo "[ERR] Expected .dict not found at {output}" >> {log}
          exit 1
        fi
        """

rule ref_bwa_index:
    input: REF
    output:
        REF + ".amb",
        REF + ".ann",
        REF + ".bwt",
        REF + ".pac",
        REF + ".sa"
    conda: f"{ENVS}/bwa_samtools.yaml"
    log:   f"{LOGS}/ref/bwa_index.log"
    threads: 1
    resources: mem_mb=2000
    shell: "mkdir -p {LOGS}/ref; bwa index {input} 2> {log}"
