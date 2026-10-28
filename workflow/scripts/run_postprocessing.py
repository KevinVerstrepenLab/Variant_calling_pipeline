#!/usr/bin/env python3
"""Execute one Snakemake analysis branch using settings from the unified config."""
import argparse
import json
from pathlib import Path

import postprocess_vcf as processing
from postprocess_vcf import Checkpointer, run_lea, run_pca_tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("pca_tree", "lea"), required=True)
    parser.add_argument("--phase", choices=("prepare", "model", "plot"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--settings", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--mem-mb", type=int, default=16000)
    args = parser.parse_args()
    settings = json.loads(args.settings)
    settings["_phase"] = args.phase
    processing.PLINK_RESOURCES = ["--memory", str(args.mem_mb), "--threads", str(args.threads)]
    args.outdir.mkdir(parents=True, exist_ok=True)
    # Native model/plot jobs remove their own products before executing. Force
    # preparation once, rather than repeating it in every downstream phase.
    checkpoint = Checkpointer(args.outdir / ".checkpoints",
                              force=settings.get("force", False) and args.phase == "prepare")
    function = run_pca_tree if args.stage == "pca_tree" else run_lea
    function(args.input.resolve(), args.outdir.resolve(), args.prefix, args.threads,
             settings, Path(__file__).resolve().parent, Path.cwd(),
             args.outdir / "commands.log", checkpoint)


if __name__ == "__main__":
    main()
