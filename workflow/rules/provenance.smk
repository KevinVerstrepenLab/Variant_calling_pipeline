rule reference_fingerprint:
    input:
        ref=REF
    output:
        sha256=f"{OUT}/provenance/reference.sha256"
    conda: str(ENVS / "gatk.yaml")
    log:
        f"{LOGS}/provenance/reference_fingerprint.log"
    threads: 1
    resources:
        mem_mb=256
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{output.sha256}')" "$(dirname '{log}')"
        sha256sum '{input.ref}' | awk '{{print $1}}' > '{output.sha256}'
        printf 'reference=%s\nsha256=%s\n' '{input.ref}' "$(cat '{output.sha256}')" > '{log}'
        """


rule gvcf_provenance:
    input:
        gvcfs=[GVCF_PATHS[sample] for sample in GVCF_SAMPLES],
        indexes=[GVCF_INDEXES[sample] for sample in GVCF_SAMPLES],
        ref=REF,
        fai=rules.ref_fai.output[0],
        reference_sha256=rules.reference_fingerprint.output.sha256,
        validator=str(WORKFLOW / "scripts" / "validate_gvcf_provenance.py")
    output:
        manifest=f"{OUT}/provenance/cohort.validated.json"
    conda: str(ENVS / "gatk.yaml")
    log:
        f"{LOGS}/provenance/cohort.log"
    params:
        samples=GVCF_SAMPLES,
        metadata=expand(f"{PROVENANCE}/{{sample}}.meta.json", sample=GVCF_SAMPLES),
        previous_metadata=PREVIOUS_GVCF_METADATA,
        logs_dir=lambda w: f"{LOGS}/provenance",
        allow_header_only=str(
            config.get("gvcf_reuse", {}).get("allow_legacy_header_match", False)
        ).lower()
    threads: config.get("gvcf_reuse", {}).get("validation_threads", 8)
    resources:
        mem_mb=lambda w, threads: threads * 512
    script:
        "../scripts/validate_cohort_provenance.py"
