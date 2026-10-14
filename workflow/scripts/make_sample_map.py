"""Assemble the validated cohort for GenomicsDBImport."""
from pathlib import Path
import json

Path(snakemake.output.map).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
expected_reference_sha256 = Path(snakemake.input.reference_sha256).read_text().strip()
manifest = json.loads(Path(snakemake.input.manifest).read_text())
if manifest.get("reference_sha256") != expected_reference_sha256:
    raise ValueError("Cohort manifest reference does not match the configured reference")
metadata_by_path = {entry["gvcf"]: entry for entry in manifest["samples"]}
expected_paths = {str(Path(path).resolve()) for path in snakemake.input.gvcfs}
if set(metadata_by_path) != expected_paths or len(metadata_by_path) != len(manifest["samples"]):
    raise ValueError("Cohort manifest must identify every selected GVCF exactly once")
rows = []
seen_samples = set()
for gvcf in snakemake.input.gvcfs:
    metadata = metadata_by_path[str(Path(gvcf).resolve())]
    sample = metadata["sample"]
    if sample in seen_samples:
        raise ValueError(f"Duplicate cohort sample: {sample}")
    if metadata.get("reference_sha256") != expected_reference_sha256:
        raise ValueError(f"Reference mismatch for reusable GVCF {gvcf}")
    seen_samples.add(sample)
    rows.append((sample, gvcf))
with open(snakemake.output.map, "w") as handle:
    for sample, gvcf in rows:
        handle.write(f"{sample}\t{gvcf}\n")
with open(snakemake.log[0], "w") as handle:
    handle.write(f"samples={len(rows)}\n")
    for sample, gvcf in rows:
        handle.write(f"{sample}\t{gvcf}\n")
