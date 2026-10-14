"""Validate a cohort in one process, with bounded parallel file readers."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace


def validate_cohort(job):
    spec = importlib.util.spec_from_file_location('gvcf_validator', job.input.validator)
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    log_path = Path(job.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    samples = list(job.params.samples)
    gvcfs = list(job.input.gvcfs)
    outputs = list(job.params.metadata)
    if not len(samples) == len(gvcfs) == len(outputs):
        raise ValueError('Cohort samples, GVCFs and metadata outputs must correspond')

    def validate_one(sample, gvcf, output):
        per_sample = SimpleNamespace(
            input=SimpleNamespace(
                gvcf=gvcf, ref=job.input.ref, fai=job.input.fai,
                reference_sha256=job.input.reference_sha256),
            output=SimpleNamespace(metadata=output),
            params=SimpleNamespace(
                allow_header_only=job.params.allow_header_only,
                previous_metadata=job.params.previous_metadata.get(sample, {})),
            wildcards=SimpleNamespace(sample=sample),
            log=[str(Path(job.params.logs_dir) / f'{sample}.log')],
        )
        validator.validate(per_sample)

    with log_path.open('w', buffering=1) as log:
        log.write(f'Validating {len(samples)} samples with {job.threads} readers\n')
        with ThreadPoolExecutor(max_workers=max(1, job.threads)) as pool:
            pending = {
                pool.submit(validate_one, sample, gvcf, output): sample
                for sample, gvcf, output in zip(samples, gvcfs, outputs)
            }
            completed = 0
            for future in as_completed(pending):
                sample = pending[future]
                try:
                    future.result()
                except Exception as error:
                    for other in pending:
                        other.cancel()
                    log.write(f'FAILED {sample}: {error}\n')
                    raise ValueError(f'Provenance validation failed for {sample}: {error}') from error
                completed += 1
                if completed % 50 == 0 or completed == len(samples):
                    log.write(f'Validated {completed}/{len(samples)} GVCFs\n')

    manifest = {
        "schema_version": 1,
        "reference_sha256": Path(job.input.reference_sha256).read_text().strip(),
        "samples": [json.loads(Path(path).read_text()) for path in outputs],
    }
    output = Path(job.output.manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary.replace(output)


if 'snakemake' in globals():
    validate_cohort(snakemake)
