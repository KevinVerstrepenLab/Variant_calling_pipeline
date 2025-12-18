rule genomicsdbimport:
    input:
        sample_map = rules.make_sample_map.output.map
    output:
        db = directory(f"{GENOMICSDB}/{{chunk}}")   # chunk is sanitized now
    params:
        interval = lambda w: INTERVAL_BY_TAG[w.chunk],
        reader_threads = config.get("joint", {}).get("reader_threads", 4),
        tmpdir = TMP_DIR
    conda: f"{ENVS}/gatk.yaml"
    log: f"{LOGS}/genomicsdb/{{chunk}}.log"
    threads: 1
    resources:
        mem_mb = config.get("genomicsdb_mem_mb", 16000)
    shell:
        r"""
        mkdir -p {GENOMICSDB}
        gatk GenomicsDBImport \
          --sample-name-map '{input.sample_map}' \
          --genomicsdb-workspace-path '{output.db}' \
          --overwrite-existing-genomicsdb-workspace true \
          -L '{params.interval}' \
          --reader-threads {params.reader_threads} \
          --tmp-dir {params.tmpdir} 2> {log}
        """
