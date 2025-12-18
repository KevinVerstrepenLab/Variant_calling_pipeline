# fastp adapter+quality trimming with auto adapter detection (good for BGI DNBseq)
rule fastp:
    input:
        R1=lambda w: PAIRS[w.sample][0],
        R2=lambda w: PAIRS[w.sample][1]
    output:
        r1 = f"{TRIM}/{{sample}}_R1.paired.fq.gz",
        r2 = f"{TRIM}/{{sample}}_R2.paired.fq.gz",
        html = f"{TRIM}/{{sample}}.fastp.html",
        json = f"{TRIM}/{{sample}}.fastp.json"
    conda: f"{ENVS}/fastp.yaml"
    log: f"{LOGS}/fastp/{{sample}}.log"
    threads: 16
    resources:
        mem_mb = config.get("fastp_mem_mb", 4000)
    log:
        f"logs/fastp/{{sample}}.log"
    params:
        min_len = config.get("fastp_min_len", 36),
        cut_mean_q = config.get("fastp_cut_mean_quality", 15),
        # Build optional flags here (strings), no Python in shell:
        flags = lambda w: " ".join([
            "--cut_front" if config.get("fastp_cut_front", True) else "",
            "--cut_tail"  if config.get("fastp_cut_tail",  True) else "",
        ]).strip()
    shell:
        r"""
        mkdir -p {TRIM} {LOGS}/fastp
        fastp \
          -i {input.R1} -I {input.R2} \
          -o {output.r1} -O {output.r2} \
          --detect_adapter_for_pe \
          --thread {threads} \
          --length_required {params.min_len} \
          {params.flags} \
          --cut_mean_quality {params.cut_mean_q} \
          --html {output.html} \
          --json {output.json} \
          2> {log}
        """
