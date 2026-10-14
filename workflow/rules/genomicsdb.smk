rule genomicsdbimport:
    input:
        sample_map=rules.make_sample_map.output.map
    output:
        db=directory(f"{GENOMICSDB}/{{chunk}}")
    conda: str(ENVS / "gatk.yaml")
    threads: config.get("joint", {}).get("reader_threads", 4)
    resources:
        mem_mb=config.get("genomicsdb_mem_mb", 16000)
    log:
        f"{LOGS}/genomicsdb/{{chunk}}.log"
    params:
        genomicsdb_dir=lambda w, output: str(Path(output[0]).parent),
        batch_size=config.get("joint", {}).get("batch_size", 50),
        consolidate="true" if config.get("joint", {}).get("consolidate", True) else "false",
        interval=lambda w: INTERVAL_BY_TAG[w.chunk],
        tmpdir=TMP_DIR,
        java_opts=lambda w, resources: f"-Xms1g -Xmx{max(1, int(resources.mem_mb * 0.8 / 1024))}g"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.genomicsdb_dir} "$(dirname '{log}')" '{params.tmpdir}'
        # A changed sample map means the interval database must be rebuilt.
        rm -rf '{output.db}'
        gatk --java-options '{params.java_opts}' GenomicsDBImport           --sample-name-map '{input.sample_map}'           --genomicsdb-workspace-path '{output.db}'           -L '{params.interval}'           --reader-threads {threads} --batch-size {params.batch_size} --consolidate {params.consolidate} --validate-sample-name-map true           --tmp-dir '{params.tmpdir}'           > '{log}' 2>&1
        """
