# Build a GATK sample map: "<sample>\t<path_to_gvcf>"
from pathlib import Path


GVCF_SUFFIX = ".g.vcf.gz"


def sample_from_gvcf(path):
    name = Path(path).name
    if not name.endswith(GVCF_SUFFIX):
        raise ValueError(f"Unexpected GVCF filename: {path}")
    return name[:-len(GVCF_SUFFIX)]


def discovered_gvcf_samples():
    """
    Include:
      1. samples from the current run, i.e. SAMPLES
      2. any existing *.g.vcf.gz files already present in GVCF
    """
    existing = [
        sample_from_gvcf(p)
        for p in Path(GVCF).glob(f"*{GVCF_SUFFIX}")
    ]

    return sorted(set(SAMPLES).union(existing))


def all_gvcfs(wildcards=None):
    return [
        f"{GVCF}/{sample}{GVCF_SUFFIX}"
        for sample in discovered_gvcf_samples()
    ]


def all_gvcf_indexes(wildcards=None):
    return [
        f"{gvcf}.tbi"
        for gvcf in all_gvcfs()
    ]


rule make_sample_map:
    input:
        gvcfs=all_gvcfs,
        idxs=all_gvcf_indexes
    output:
        map=f"{JOINT}/samples.map"
    run:
        import os

        os.makedirs(JOINT, exist_ok=True)

        with open(output.map, "w") as fh:
            for gvcf in sorted(input.gvcfs):
                sample = sample_from_gvcf(gvcf)
                fh.write(f"{sample}\t{gvcf}\n")