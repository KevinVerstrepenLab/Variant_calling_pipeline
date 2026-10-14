"""Shared workflow path, interval and filtering helpers."""
from pathlib import Path
import os
import re
import gzip
from snakemake.exceptions import WorkflowError

GVCF_SUFFIX = ".g.vcf.gz"

def sanitize(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", value)


def _load_contigs_from_dict(dict_path):
    contigs = []
    with open(dict_path) as handle:
        for line in handle:
            if line.startswith("@SQ"):
                parts = dict(x.split(":", 1) for x in line.strip().split("\t")[1:])
                contigs.append((parts["SN"], int(parts["LN"])))
    return contigs


def _load_contigs_from_fai(fai_path):
    contigs = []
    with open(fai_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields:
                contigs.append((fields[0], int(fields[1])))
    return contigs


def _scan_fasta_lengths(fasta_path):
    """Last-resort parse-time reader; supports gzip-compressed FASTA."""
    contigs = []
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    name, length = None, 0
    with opener(fasta_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        contigs.append((name, length))
    if not contigs:
        raise WorkflowError(f"Failed to derive contigs from FASTA: {fasta_path}")
    return contigs


def _expanduser(path):
    return Path(os.path.expanduser(str(path))).resolve()


def trimmed_reads(sample, LAYOUT, TRIM):
    if LAYOUT[sample] == "paired":
        return [
            f"{TRIM}/{sample}_R1.paired.fq.gz",
            f"{TRIM}/{sample}_R2.paired.fq.gz",
        ]
    return [f"{TRIM}/{sample}.single.fq.gz"]


def sample_from_gvcf(path):
    name = Path(path).name
    if not name.endswith(GVCF_SUFFIX):
        raise ValueError(f"Unexpected GVCF filename: {path}")
    return name[:-len(GVCF_SUFFIX)]

def ad_expr(config):
    ad = config.get("ad_filter", {})
    min_dp = int(ad.get("min_dp", 5))
    lo = float(ad.get("ambig_low", 0.1))
    hi = float(ad.get("ambig_high", 0.9))
    # [:allele] selects that allele across samples; single | and & retain
    # sample-specific masks instead of selecting every genotype at the site.
    return (
        f'FMT/DP<{min_dp} | '
        f'((FMT/AD[:0]+FMT/AD[:1])>0 & '
        f'FMT/AD[:1]/(FMT/AD[:0]+FMT/AD[:1])>{lo} & '
        f'FMT/AD[:1]/(FMT/AD[:0]+FMT/AD[:1])<{hi})'
    )
