# Separate rules avoid untracked placeholder outputs while allowing mixed layouts.
rule fastp_paired:
    input:
        r1=lambda w: READS[w.sample]["R1"],
        r2=lambda w: READS[w.sample]["R2"]
    output:
        r1=f"{TRIM}/{{sample}}_R1.paired.fq.gz",
        r2=f"{TRIM}/{{sample}}_R2.paired.fq.gz",
        html=f"{TRIM}/{{sample}}.fastp.html",
        json=f"{TRIM}/{{sample}}.fastp.json"
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in PAIRED_SAMPLES) if PAIRED_SAMPLES else r"(?!)"
    conda: str(ENVS / "fastp.yaml")
    threads: config.get("fastp_threads", 16)
    resources:
        mem_mb=config.get("fastp_mem_mb", 4000)
    log:
        f"{LOGS}/fastp/{{sample}}.log"
    params:
        trim_dir=lambda w, output: str(Path(output[0]).parent),
        min_len=config.get("fastp_min_len", 36),
        cut_mean_q=config.get("fastp_cut_mean_quality", 15),
        flags=lambda w: " ".join(flag for flag in [
            "--cut_front" if config.get("fastp_cut_front", True) else "",
            "--cut_tail" if config.get("fastp_cut_tail", True) else "",
        ] if flag)
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.trim_dir} "$(dirname '{log}')"
        fastp           -i '{input.r1}' -I '{input.r2}'           -o '{output.r1}' -O '{output.r2}'           --detect_adapter_for_pe           --thread {threads}           --length_required {params.min_len}           {params.flags}           --cut_mean_quality {params.cut_mean_q}           --html '{output.html}'           --json '{output.json}'           > '{log}' 2>&1
        """


rule fastp_single:
    input:
        r1=lambda w: READS[w.sample]["R1"]
    output:
        r1=f"{TRIM}/{{sample}}.single.fq.gz",
        html=f"{TRIM}/{{sample}}.fastp.html",
        json=f"{TRIM}/{{sample}}.fastp.json"
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SINGLE_SAMPLES) if SINGLE_SAMPLES else r"(?!)"
    conda: str(ENVS / "fastp.yaml")
    threads: config.get("fastp_threads", 16)
    resources:
        mem_mb=config.get("fastp_mem_mb", 4000)
    log:
        f"{LOGS}/fastp/{{sample}}.log"
    params:
        trim_dir=lambda w, output: str(Path(output[0]).parent),
        min_len=config.get("fastp_min_len", 36),
        cut_mean_q=config.get("fastp_cut_mean_quality", 15),
        flags=lambda w: " ".join(flag for flag in [
            "--cut_front" if config.get("fastp_cut_front", True) else "",
            "--cut_tail" if config.get("fastp_cut_tail", True) else "",
        ] if flag)
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.trim_dir} "$(dirname '{log}')"
        fastp           -i '{input.r1}'           -o '{output.r1}'           --thread {threads}           --length_required {params.min_len}           {params.flags}           --cut_mean_quality {params.cut_mean_q}           --html '{output.html}'           --json '{output.json}'           > '{log}' 2>&1
        """
