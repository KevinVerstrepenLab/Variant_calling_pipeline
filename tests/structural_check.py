#!/usr/bin/env python3
"""Dependency-free structural checks for the public workflow package."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
errors: list[str] = []

required = [
    "README.md",
    "VERSION",
    "config/config.yaml",
    "config/samples_examples.tsv",
    "tests/minimal/config.yaml",
    "tests/minimal/samples.tsv",
    "workflow/configuration.py",
    "workflow/rules/postprocessing.smk",
    "workflow/rules/coverage.smk",
    "docs/RELEASE_NOTES_v1.0.0.md",
    "docs/UPDATE_v1.0.0.md",
    "CHANGELOG.md",
    "RELEASE_CHECKLIST.md",
    "STRUCTURAL_REVIEW.md",
    "workflow/Snakefile",
    "workflow/rules/merge.smk",
    "workflow/rules/provenance.smk",
    "workflow/scripts/validate_gvcf_provenance.py",
    "docs/gvcf_reuse.md",
    ".github/workflows/snakemake-ci.yml",
]
for relative in required:
    if not (ROOT / relative).is_file():
        errors.append(f"missing required file: {relative}")

# Inspect distributed source only, excluding runtime caches and local datasets.
source_files = [p for directory in ("workflow", "docs", "tests", ".github")
                for p in (ROOT / directory).rglob("*")
                if p.is_file() and "__pycache__" not in p.parts
                and "postprocessing_bundle" not in p.parts]
source_files += list(ROOT.glob("*.md"))
source_files += [ROOT / "config/config.yaml", ROOT / "config/samples_examples.tsv"]
private_path = re.compile(r"/data\d+/u\d+|/home/[^/]+")
for path in source_files:
    if path.name.endswith("_old") or path.suffix == ".unix":
        continue  # preserved local backups are not part of the distributed source
    if path.suffix not in {".smk", ".py", ".R", ".yaml", ".yml", ".md", ".tsv"}:
        continue
    if path == Path(__file__).resolve():
        continue
    if private_path.search(path.read_text(errors="replace")):
        errors.append(f"private absolute path found: {path.relative_to(ROOT)}")

for path in (ROOT / "workflow/rules").glob("*.smk"):
    text = path.read_text()
    # Catch repeated top-level rule directives within a single rule block.
    blocks = re.split(r"(?m)^[ \t]*rule\s+", text)[1:]
    for block in blocks:
        rule_name = block.split(":", 1)[0].strip()
        for directive in ("input", "output", "params", "threads", "resources", "log", "shell", "script", "run", "conda"):
            count = len(re.findall(rf"(?m)^[ \t]+{directive}:\s*", block))
            if count > 1:
                errors.append(f"duplicate {directive}: in {path.name}:{rule_name}")
        if re.search(r"(?m)^[ \t]+shell:\s*", block) and "set -euo pipefail" not in block:
            errors.append(f"shell rule lacks strict mode: {path.name}:{rule_name}")

for path in (ROOT / "workflow/rules").glob("*.smk"):
    for script in re.findall(r'script:\s*"([^"]+)"', path.read_text()):
        if not (path.parent / script).is_file():
            errors.append(f"broken script path: {path.name}: {script}")

snakefile = (ROOT / "workflow/Snakefile").read_text()
for include in re.findall(r'include:\s*"([^"]+)"', snakefile):
    if not (ROOT / "workflow" / include).is_file():
        errors.append(f"broken include: {include}")

if errors:
    print("STRUCTURAL CHECK FAILED")
    for error in errors:
        print(f"- {error}")
    sys.exit(1)
print("STRUCTURAL CHECK PASSED")
