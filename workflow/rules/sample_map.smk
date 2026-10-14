# Build a GATK sample map from current samples and complete reusable GVCFs.
from pathlib import Path

GVCF_SUFFIX = ".g.vcf.gz"



rule make_sample_map:
    input:
        gvcfs=[GVCF_PATHS[sample] for sample in GVCF_SAMPLES],
        idxs=[GVCF_INDEXES[sample] for sample in GVCF_SAMPLES],
        manifest=rules.gvcf_provenance.output.manifest,
        reference_sha256=rules.reference_fingerprint.output.sha256
    output:
        map=f"{JOINT}/samples.map"
    log:
        f"{LOGS}/joint/sample_map.log"
    conda: str(ENVS / "gatk.yaml")
    script:
        "../scripts/make_sample_map.py"
