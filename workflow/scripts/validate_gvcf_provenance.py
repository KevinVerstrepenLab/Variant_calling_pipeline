#!/usr/bin/env python3
"""Validate that a reusable GVCF was produced against the configured reference."""

import gzip
import hashlib
import json
import re
import shlex
from pathlib import Path
from urllib.parse import unquote, urlparse


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fai(path: Path) -> list[tuple[str, int]]:
    contigs = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                contigs.append((fields[0], int(fields[1])))
    return contigs


def parse_contig(line: str) -> tuple[str, int] | None:
    if not line.startswith("##contig=<"):
        return None
    body = line[len("##contig=<") :].rstrip("\n>")
    fields = {}
    for item in body.split(","):
        if "=" in item:
            key, value = item.split("=", 1)
            fields[key] = value.strip('"')
    if "ID" not in fields or "length" not in fields:
        return None
    return fields["ID"], int(fields["length"])


def read_gvcf_header(path: Path) -> tuple[str | None, list[tuple[str, int]]]:
    reference_uri = None
    contigs = []
    samples = None
    command_reference = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##reference="):
                reference_uri = line.split("=", 1)[1].strip()
            if line.startswith("##GATKCommandLine=<ID=HaplotypeCaller,"):
                match = re.search(r'CommandLine="(.*?)",Version=', line)
                if match:
                    arguments = shlex.split(match.group(1))
                    for index, argument in enumerate(arguments[:-1]):
                        if argument in {"--reference", "-R"}:
                            command_reference = arguments[index + 1]
            parsed = parse_contig(line)
            if parsed:
                contigs.append(parsed)
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                # Consume the stream to detect truncated gzip files, including CRC failures.
                if str(path).endswith(".gz"):
                    for _ in handle:
                        pass
                break
    if samples is None:
        raise ValueError(f"Missing #CHROM header in {path}")
    return reference_uri or command_reference, contigs, samples


def local_reference_from_uri(uri: str | None) -> Path | None:
    if not uri:
        return None
    if uri.startswith("file:"):
        parsed = urlparse(uri)
        return Path(unquote(parsed.path))
    candidate = Path(uri)
    return candidate if candidate.is_absolute() else None


def validate(snakemake):
    gvcf = Path(snakemake.input.gvcf)
    reference = Path(snakemake.input.ref)
    fai = Path(snakemake.input.fai)
    reference_sha256 = Path(snakemake.input.reference_sha256).read_text().strip()
    metadata_path = Path(snakemake.output.metadata)
    log_path = Path(snakemake.log[0])
    allow_header_only = str(snakemake.params.allow_header_only).lower() == "true"

    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    reference_uri, gvcf_contigs, samples = read_gvcf_header(gvcf)
    if samples != [snakemake.wildcards.sample]:
        raise ValueError(f"GVCF sample mismatch: expected {snakemake.wildcards.sample!r}, found {samples}")
    gvcf_sha256 = sha256_file(gvcf)
    # Snakemake removes old outputs before rerunning a job. The workflow captures
    # the previous sidecar before that removal so legacy adoption survives reruns.
    previous = getattr(snakemake.params, "previous_metadata", None)
    if previous is None:
        previous = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    previous_valid = (
        previous.get("sample") == snakemake.wildcards.sample
        and previous.get("reference_sha256") == reference_sha256
        and previous.get("gvcf_sha256") == gvcf_sha256
    )
    expected_contigs = read_fai(fai)
    if gvcf_contigs != expected_contigs:
        raise ValueError(
            f"GVCF/reference contig mismatch for {gvcf}. "
            "The GVCF cannot be reused with the configured reference."
        )

    validation_mode = None
    header_reference = local_reference_from_uri(reference_uri)
    header_reference_sha256 = None
    if header_reference and header_reference.is_file():
        header_reference_sha256 = sha256_file(header_reference)
        if header_reference_sha256 != reference_sha256:
            raise ValueError(
                f"Reference SHA-256 mismatch for reusable GVCF {gvcf}. "
                f"GVCF header reference: {header_reference}; configured reference: {reference}."
            )
        validation_mode = "header_reference_sha256"
    elif previous_valid:
        validation_mode = "validated_sidecar"
    elif allow_header_only:
        validation_mode = "contig_header_only"
    else:
        raise ValueError(
            f"Cannot prove that {gvcf} used the configured reference. Its ##reference path "
            "is absent or inaccessible. To adopt legacy GVCFs after independently confirming "
            "the reference, set gvcf_reuse.allow_legacy_header_match: true once."
        )

    metadata = {
        "schema_version": 2,
        "gvcf_sha256": gvcf_sha256,
        "sample": snakemake.wildcards.sample,
        "gvcf": str(gvcf.resolve()),
        "gvcf_size_bytes": gvcf.stat().st_size,
        "reference": str(reference.resolve()),
        "reference_sha256": reference_sha256,
        "gvcf_header_reference": reference_uri,
        "gvcf_header_reference_sha256": header_reference_sha256,
        "validation_mode": validation_mode,
        "contigs": [{"name": name, "length": length} for name, length in expected_contigs],
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    log_path.write_text(
        f"Validated {gvcf}\nreference_sha256={reference_sha256}\n"
        f"validation_mode={validation_mode}\n"
    )


if "snakemake" in globals():
    validate(snakemake)
