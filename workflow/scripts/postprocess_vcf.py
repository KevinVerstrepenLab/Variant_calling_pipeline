#!/usr/bin/env python3
"""Post-processing implementation shared by the two Snakemake analysis branches.

Settings and input paths are supplied by run_postprocessing.py from the unified
Snakemake configuration. Internal checkpoints retain expensive intermediates
when only downstream plotting options change.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

CHECKPOINT_VERSION = 11
SMALL_FILE_HASH_LIMIT = 25 * 1024 * 1024
PLINK_RESOURCES: List[str] = []


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def as_bool(x: Any, default: bool = False) -> bool:
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    return s in {"1", "true", "yes", "y", "on"}


def q(path_or_arg: Any) -> str:
    return shlex.quote(str(path_or_arg))


def relpath(value: Optional[Any], base: Path) -> Optional[Path]:
    if value is None or str(value).strip() == "":
        return None
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else base / path


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def require_tools(tools: Iterable[str]) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        die("missing required executables in PATH: " + ", ".join(missing))


def stable_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def path_fingerprint(path: Optional[Path]) -> Any:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    stat = p.stat()
    info: Dict[str, Any] = {
        "path": str(p.resolve()),
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "type": "dir" if p.is_dir() else "file",
    }
    if p.is_file() and stat.st_size <= SMALL_FILE_HASH_LIMIT:
        info["sha256"] = sha256_file(p)
    return info


def normalize_for_signature(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value.resolve()) if value.exists() else str(value)
    if isinstance(value, Mapping):
        return {str(k): normalize_for_signature(v) for k, v in sorted(value.items(), key=lambda x: str(x[0]))}
    if isinstance(value, (list, tuple)):
        return [normalize_for_signature(v) for v in value]
    return value


def output_ready(path: Path) -> bool:
    p = Path(path)
    if p.is_dir():
        return p.exists()
    if not p.exists():
        return False
    # Some checkpointed outputs are valid even when empty.
    # Example: <prefix>.remove_samples.txt is intentionally empty when no samples
    # exceed the missingness threshold. Requiring size > 0 made such successful
    # filter steps fail after vcftools completed correctly.
    if p.name.endswith(".remove_samples.txt"):
        return True
    return p.stat().st_size > 0


def outputs_ready(outputs: Sequence[Path]) -> bool:
    return all(output_ready(Path(p)) for p in outputs)


def run(cmd: Sequence[str], log: Optional[Path] = None, cwd: Optional[Path] = None) -> None:
    """Run command, stream stdout/stderr to console and optional log file."""
    cmd = list(cmd)
    if cmd and cmd[0] in {"plink", "plink2"}:
        for flag, value in zip(PLINK_RESOURCES[::2], PLINK_RESOURCES[1::2]):
            if flag not in cmd:
                cmd.extend([flag, value])
    printable = " ".join(q(x) for x in cmd)
    print(f"\n[command] {printable}")
    log_handle = None
    try:
        if log is not None:
            mkdir(log.parent)
            log_handle = open(log, "a", encoding="utf-8")
            log_handle.write(f"\n# {dt.datetime.now().isoformat(timespec='seconds')}\n{printable}\n")
            log_handle.flush()
        proc = subprocess.Popen(
            list(map(str, cmd)),
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            if log_handle:
                log_handle.write(line)
        rc = proc.wait()
        if rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)
    finally:
        if log_handle:
            log_handle.close()


def index_vcf(vcf: Path, threads: int, log: Optional[Path]) -> None:
    run(["bcftools", "index", "-f", "--threads", str(threads), str(vcf)], log=log)


class Checkpointer:
    def __init__(self, checkpoint_dir: Path, enabled: bool = True, force: bool = False) -> None:
        self.checkpoint_dir = checkpoint_dir
        self.enabled = enabled
        self.force = force
        mkdir(checkpoint_dir)

    @staticmethod
    def safe_name(name: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")

    def stamp_path(self, name: str) -> Path:
        return self.checkpoint_dir / f"{self.safe_name(name)}.json"

    def payload(self, name: str, inputs: Mapping[str, Optional[Path]], outputs: Sequence[Path], params: Mapping[str, Any]) -> Dict[str, Any]:
        return {
            "checkpoint_version": CHECKPOINT_VERSION,
            "step": name,
            "inputs": {k: path_fingerprint(v) for k, v in sorted(inputs.items())},
            "outputs": [str(Path(p).resolve()) for p in outputs],
            "params": normalize_for_signature(params),
        }

    def is_current(self, name: str, inputs: Mapping[str, Optional[Path]], outputs: Sequence[Path], params: Mapping[str, Any]) -> bool:
        if not self.enabled or self.force:
            return False
        if not outputs_ready(outputs):
            return False
        stamp = self.stamp_path(name)
        if not stamp.exists():
            return False
        payload = self.payload(name, inputs, outputs, params)
        expected = sha256_text(stable_json(payload))
        try:
            with open(stamp, "r", encoding="utf-8") as handle:
                saved = json.load(handle)
        except Exception:
            return False
        return saved.get("hash") == expected

    def write_stamp(self, name: str, inputs: Mapping[str, Optional[Path]], outputs: Sequence[Path], params: Mapping[str, Any]) -> None:
        payload = self.payload(name, inputs, outputs, params)
        stamp = self.stamp_path(name)
        record = {
            "hash": sha256_text(stable_json(payload)),
            "created_at": dt.datetime.now().isoformat(timespec="seconds"),
            "payload": payload,
        }
        with open(stamp, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.write("\n")

    def run(self, name: str, inputs: Mapping[str, Optional[Path]], outputs: Sequence[Path], params: Mapping[str, Any], action: Any) -> None:
        outputs = [Path(p) for p in outputs]
        if self.is_current(name, inputs, outputs, params):
            print(f"\n[checkpoint] {name}: up to date; skipping")
            return
        if self.enabled:
            reason = "forced" if self.force else "missing outputs or changed input/parameters"
            print(f"\n[checkpoint] {name}: running ({reason})")
        else:
            print(f"\n[checkpoint] {name}: checkpointing disabled; running")
        action()
        missing = [str(p) for p in outputs if not output_ready(p)]
        if missing:
            die("step did not create expected output(s): " + ", ".join(missing))
        if self.enabled:
            self.write_stamp(name, inputs, outputs, params)

    def clear(self) -> None:
        if self.checkpoint_dir.exists():
            for p in self.checkpoint_dir.glob("*.json"):
                p.unlink()

    def list(self) -> None:
        print(f"Checkpoint directory: {self.checkpoint_dir}")
        if not self.checkpoint_dir.exists():
            print("No checkpoint directory found.")
            return
        stamps = sorted(self.checkpoint_dir.glob("*.json"))
        if not stamps:
            print("No checkpoints found.")
            return
        for stamp in stamps:
            try:
                with open(stamp, "r", encoding="utf-8") as handle:
                    data = json.load(handle)
                payload = data.get("payload", {})
                print(f"- {payload.get('step', stamp.stem)} | {data.get('created_at', 'unknown time')}")
            except Exception:
                print(f"- {stamp.name} | unreadable")


def maybe_add_path_arg(args: List[str], flag: str, path: Optional[Path]) -> None:
    if path is not None:
        args.extend([flag, str(path)])


def pgen_outputs(prefix_path: Path) -> List[Path]:
    return [Path(str(prefix_path) + ext) for ext in [".pgen", ".pvar", ".psam"]]


def bed_outputs(prefix_path: Path) -> List[Path]:
    return [Path(str(prefix_path) + ext) for ext in [".bed", ".bim", ".fam"]]


def run_pca_tree(input_vcf: Path, outdir: Path, prefix: str, threads: int, cfg: Mapping[str, Any], scripts_dir: Path, config_base: Path, log: Path, ckpt: Checkpointer) -> None:
    enabled = as_bool(cfg.get("enabled"), default=False)
    if not enabled:
        print("\n[pca_tree] disabled")
        return

    require_tools(["plink2", "plink", "Rscript"])
    mkdir(outdir)
    pfx = outdir / prefix

    maf = str(cfg.get("maf", 0.01))
    geno = str(cfg.get("geno", 0.10))
    method = str(cfg.get("method", "nj"))
    components = int(cfg.get("pca_components", 20))
    do_pca = as_bool(cfg.get("pca"), default=True)
    do_tree = as_bool(cfg.get("tree"), default=True)

    keep = relpath(cfg.get("keep"), config_base)
    pca_sample_list = relpath(cfg.get("pca_sample_list"), config_base)
    tree_sample_list = relpath(cfg.get("tree_sample_list"), config_base)
    rename = relpath(cfg.get("rename"), config_base)
    clades = relpath(cfg.get("clades"), config_base)
    if clades is not None and (not clades.exists() or clades.stat().st_size == 0):
        print(f"[pca_tree] clades file not found or empty; continuing without clade colouring: {clades}")
        clades = None
    clade_colors = relpath(cfg.get("clade_colors", cfg.get("clade_palette", cfg.get("palette"))), config_base)
    if clade_colors is not None and (not clade_colors.exists() or clade_colors.stat().st_size == 0):
        print(f"[pca_tree] clade_colors file not found or empty; generating palette automatically: {clade_colors}")
        clade_colors = None
    outgroup = cfg.get("outgroup", None)
    layout = str(cfg.get("layout", "fan")).strip()
    tree_plot_prefix = outdir / f"{prefix}.tree-{layout}"
    branch_length = as_bool(cfg.get("branch_length", cfg.get("branch.length", True)), default=True)
    width = str(cfg.get("width", 18))
    height = str(cfg.get("height", 18))
    dpi = str(cfg.get("dpi", 300))
    label_tips = as_bool(cfg.get("label_tips"), default=False)
    label_size = str(cfg.get("label_size", cfg.get("label.size", 2.0)))
    pca_label_points = as_bool(cfg.get("pca_label_points"), default=label_tips)
    pca_label_size = str(cfg.get("pca_label_size", label_size))
    tip_size = str(cfg.get("tip_size", cfg.get("tip.size", cfg.get("tips_size", cfg.get("tips.size", 1.2)))))
    color_branches = as_bool(cfg.get("color_branches", cfg.get("colour_branches", cfg.get("branch_color", cfg.get("branch_colour", True)))), default=True)
    color_tip_labels = as_bool(cfg.get("color_tip_labels", cfg.get("colour_tip_labels", cfg.get("color_labels", cfg.get("colour_labels", True)))), default=True)
    snps_only = as_bool(cfg.get("snps_only"), default=True)
    biallelic = as_bool(cfg.get("biallelic"), default=True)

    def action_import() -> None:
        args = [
            "plink2", "--vcf", str(input_vcf), "--allow-extra-chr", "--threads", str(threads),
            "--maf", maf, "--geno", geno,
        ]
        if snps_only:
            args.append("--snps-only")
        if biallelic:
            args.extend(["--max-alleles", "2", "--min-alleles", "2"])
        maybe_add_path_arg(args, "--keep", keep)
        args.extend(["--make-pgen", "--out", str(pfx)])
        run(args, log=log)

    ckpt.run(
        f"pca_tree.import.{prefix}",
        inputs={"input_vcf": input_vcf, "keep": keep},
        outputs=pgen_outputs(pfx),
        params={"threads": threads, "maf": maf, "geno": geno, "snps_only": snps_only, "biallelic": biallelic},
        action=action_import,
    )

    uniq = outdir / f"{prefix}.uniq"

    def action_unique() -> None:
        run([
            "plink2", "--pfile", str(pfx), "--allow-extra-chr", "--threads", str(threads),
            "--set-all-var-ids", "@:#:$r:$a", "--rm-dup", "force-first", "--make-pgen", "--out", str(uniq)
        ], log=log)

    ckpt.run(
        f"pca_tree.unique_ids.{prefix}",
        inputs={"pgen": Path(str(pfx) + ".pgen"), "pvar": Path(str(pfx) + ".pvar"), "psam": Path(str(pfx) + ".psam")},
        outputs=pgen_outputs(uniq),
        params={"threads": threads, "var_ids": "@:#:$r:$a", "rm_dup": "force-first"},
        action=action_unique,
    )
    active = uniq

    ld = cfg.get("ld_prune", {}) or {}
    if as_bool(ld.get("enabled"), default=True):
        window = str(ld.get("window", 50))
        step = str(ld.get("step", 5))
        r2 = str(ld.get("r2", 0.2))
        prune_prefix = outdir / f"{prefix}.prune"
        prune_in = Path(str(prune_prefix) + ".prune.in")

        def action_prune_list() -> None:
            run([
                "plink2", "--pfile", str(active), "--allow-extra-chr", "--threads", str(threads),
                "--indep-pairwise", window, step, r2, "--out", str(prune_prefix)
            ], log=log)

        ckpt.run(
            f"pca_tree.prune_list.{prefix}",
            inputs={"active_pgen": Path(str(active) + ".pgen"), "active_pvar": Path(str(active) + ".pvar"), "active_psam": Path(str(active) + ".psam")},
            outputs=[prune_in],
            params={"threads": threads, "window": window, "step": step, "r2": r2},
            action=action_prune_list,
        )

        pruned = outdir / f"{prefix}.pruned"

        def action_extract_pruned() -> None:
            run([
                "plink2", "--pfile", str(active), "--allow-extra-chr", "--threads", str(threads),
                "--extract", str(prune_in), "--make-pgen", "--out", str(pruned)
            ], log=log)

        ckpt.run(
            f"pca_tree.extract_pruned.{prefix}",
            inputs={"active_pgen": Path(str(active) + ".pgen"), "active_pvar": Path(str(active) + ".pvar"), "active_psam": Path(str(active) + ".psam"), "prune_in": prune_in},
            outputs=pgen_outputs(pruned),
            params={"threads": threads},
            action=action_extract_pruned,
        )
        active = pruned

    eigenvec: Optional[Path] = None
    eigenval: Optional[Path] = None
    if do_pca:
        pca_prefix = outdir / f"{prefix}.pca"
        eigenvec = Path(str(pca_prefix) + ".eigenvec")
        eigenval = Path(str(pca_prefix) + ".eigenval")

        def action_pca() -> None:
            run([
                "plink2", "--pfile", str(active), "--allow-extra-chr", "--threads", str(threads),
                "--pca", str(components), "approx", "--out", str(pca_prefix)
            ], log=log)

        ckpt.run(
            f"pca_tree.pca.{prefix}",
            inputs={"active_pgen": Path(str(active) + ".pgen"), "active_pvar": Path(str(active) + ".pvar"), "active_psam": Path(str(active) + ".psam")},
            outputs=[eigenvec, eigenval],
            params={"threads": threads, "components": components},
            action=action_pca,
        )

    if do_tree:
        bedset = Path(str(active) + ".bedset")

        def action_make_bed() -> None:
            run([
                "plink2", "--pfile", str(active), "--allow-extra-chr", "--threads", str(threads),
                "--max-alleles", "2", "--min-alleles", "2", "--make-bed", "--out", str(bedset)
            ], log=log)

        ckpt.run(
            f"pca_tree.make_bed.{prefix}",
            inputs={"active_pgen": Path(str(active) + ".pgen"), "active_pvar": Path(str(active) + ".pvar"), "active_psam": Path(str(active) + ".psam")},
            outputs=bed_outputs(bedset),
            params={"threads": threads, "biallelic": True},
            action=action_make_bed,
        )

        dist_prefix = outdir / f"{prefix}.dist"
        mdist = Path(str(dist_prefix) + ".mdist")
        mdist_id = Path(str(dist_prefix) + ".mdist.id")

        def action_distance() -> None:
            run([
                "plink", "--bfile", str(bedset), "--allow-extra-chr", "--threads", str(threads),
                "--distance", "square", "1-ibs", "--out", str(dist_prefix)
            ], log=log)

        ckpt.run(
            f"pca_tree.distance.{prefix}",
            inputs={"bed": Path(str(bedset) + ".bed"), "bim": Path(str(bedset) + ".bim"), "fam": Path(str(bedset) + ".fam")},
            outputs=[mdist, mdist_id],
            params={"threads": threads, "distance": "square 1-ibs"},
            action=action_distance,
        )

    if cfg.get("_phase") == "prepare":
        return

    if do_pca or do_tree:
        tree_prefix = outdir / f"{prefix}.tree"
        rscript = scripts_dir / "build_tree_postprocess.R"
        tree_outputs = [
            Path(str(tree_prefix) + ".nwk"),
            Path(str(tree_prefix) + ".tree_info.tsv"),
            Path(str(tree_plot_prefix) + ".pdf"),
            Path(str(tree_plot_prefix) + ".png"),
        ] if do_tree else []
        if clades is not None:
            tree_outputs.append(Path(str(tree_prefix) + ".clade_colors.tsv"))
        if eigenvec is not None and eigenval is not None:
            tree_outputs += [Path(str(tree_prefix) + ext) for ext in [".pca_scores.tsv", ".pca.pdf", ".pca.png"]]

        def action_tree_plot() -> None:
            rargs = [
                "Rscript", str(rscript),
                "--method", method, "--out", str(tree_prefix), "--layout", layout,
                "--branch_length", "yes" if branch_length else "no",
                "--width", width, "--height", height, "--dpi", dpi,
                "--label_size", label_size, "--tip_size", tip_size,
                "--pca_label_points", "yes" if pca_label_points else "no",
                "--pca_label_size", pca_label_size,
            ]
            if do_tree:
                rargs.extend(["--dist", str(mdist), "--ids", str(mdist_id)])
            if color_branches:
                rargs.append("--color_branches")
            if color_tip_labels:
                rargs.append("--color_tip_labels")
            maybe_add_path_arg(rargs, "--rename", rename)
            maybe_add_path_arg(rargs, "--pca_sample_list", pca_sample_list)
            maybe_add_path_arg(rargs, "--tree_sample_list", tree_sample_list)
            maybe_add_path_arg(rargs, "--clades", clades)
            maybe_add_path_arg(rargs, "--clade_colors", clade_colors)
            if outgroup is not None and str(outgroup).strip():
                rargs.extend(["--outgroup", str(outgroup)])
            if label_tips:
                rargs.append("--label_tips")
            if eigenvec is not None and eigenval is not None:
                rargs.extend(["--pca_eigenvec", str(eigenvec), "--pca_eigenval", str(eigenval)])
            run(rargs, log=log)

        ckpt.run(
            f"pca_tree.plot_tree_pca.{prefix}",
            inputs={"mdist": mdist if do_tree else None, "mdist_id": mdist_id if do_tree else None, "rename": rename, "pca_sample_list": pca_sample_list, "tree_sample_list": tree_sample_list, "clades": clades, "clade_colors": clade_colors, "eigenvec": eigenvec, "eigenval": eigenval, "rscript": rscript, "palette_helper": scripts_dir / "clade_palette_utils.R"},
            outputs=tree_outputs,
            params={"method": method, "outgroup": outgroup, "layout": layout, "branch_length": branch_length, "width": width, "height": height, "dpi": dpi, "label_tips": label_tips, "label_size": label_size, "pca_label_points": pca_label_points, "pca_label_size": pca_label_size, "tip_size": tip_size, "color_branches": color_branches, "color_tip_labels": color_tip_labels, "clade_colors": str(clade_colors) if clade_colors else None},
            action=action_tree_plot,
        )

def run_lea(input_vcf: Path, outdir: Path, prefix: str, threads: int, cfg: Mapping[str, Any], scripts_dir: Path, config_base: Path, log: Path, ckpt: Checkpointer) -> None:
    enabled = as_bool(cfg.get("enabled"), default=False)
    if not enabled:
        print("\n[lea] disabled")
        return

    require_tools(["bcftools", "vcftools", "plink", "Rscript"])
    mkdir(outdir)
    maf = str(cfg.get("maf", 0.05))
    max_missing = str(cfg.get("max_missing", 0.9))
    ld = cfg.get("ld_prune", {}) or {}
    window = str(ld.get("window", 50))
    step = str(ld.get("step", 5))
    r2 = str(ld.get("r2", 0.2))
    rename = relpath(cfg.get("rename"), config_base)
    clades = relpath(cfg.get("clades"), config_base)
    if clades is not None and (not clades.exists() or clades.stat().st_size == 0):
        print(f"[lea] clades file not found or empty; continuing without clade colour annotations: {clades}")
        clades = None
    clade_colors = relpath(cfg.get("clade_colors", cfg.get("clade_palette", cfg.get("palette"))), config_base)
    if clade_colors is not None and (not clade_colors.exists() or clade_colors.stat().st_size == 0):
        print(f"[lea] clade_colors file not found or empty; generating palette automatically: {clade_colors}")
        clade_colors = None
    color_sample_labels_by_clade = as_bool(cfg.get("color_sample_labels_by_clade", cfg.get("color_labels_by_clade", True)), default=True)
    clade_color_strip = as_bool(cfg.get("clade_color_strip", cfg.get("color_strip", True)), default=True)

    raw_bial = outdir / f"{prefix}.biallelic.raw.vcf.gz"
    bial = outdir / f"{prefix}.biallelic.vcf.gz"
    filt_prefix = outdir / f"{prefix}.biallelic.maf{str(maf).replace('.', '')}.miss{str(max_missing).replace('.', '')}"
    filt_vcf_plain = Path(str(filt_prefix) + ".recode.vcf")
    filt_vcf = Path(str(filt_prefix) + ".vcf.gz")
    plink_prefix = outdir / f"{prefix}.plink"
    prune_prefix = outdir / f"{prefix}.prune"
    lea_bed = outdir / f"{prefix}.LEA"
    lea_final = outdir / f"{prefix}.LEA.final"
    pruned_vcf = outdir / f"{prefix}.pruned.vcf.gz"

    def action_biallelic() -> None:
        run(["bcftools", "view", "-m2", "-M2", "-v", "snps", "--threads", str(threads), str(input_vcf), "-Oz", "-o", str(raw_bial)], log=log)
        index_vcf(raw_bial, threads, log)

    ckpt.run(
        f"lea.biallelic_snps.{prefix}",
        inputs={"input_vcf": input_vcf},
        outputs=[raw_bial, Path(str(raw_bial) + ".csi")],
        params={"threads": threads, "biallelic": True, "variant_type": "snps"},
        action=action_biallelic,
    )

    def action_set_ids() -> None:
        run(["bcftools", "annotate", "--set-id", "%CHROM:%POS:%REF:%ALT", "--threads", str(threads), "-Oz", "-o", str(bial), str(raw_bial)], log=log)
        index_vcf(bial, threads, log)

    ckpt.run(
        f"lea.set_variant_ids.{prefix}",
        inputs={"raw_biallelic_vcf": raw_bial},
        outputs=[bial, Path(str(bial) + ".csi")],
        params={"threads": threads, "set_id": "%CHROM:%POS:%REF:%ALT"},
        action=action_set_ids,
    )

    def action_vcftools_filter() -> None:
        run(["vcftools", "--gzvcf", str(bial), "--maf", maf, "--max-missing", max_missing, "--recode", "--recode-INFO-all", "--out", str(filt_prefix)], log=log)
        run(["bcftools", "view", "--threads", str(threads), "-Oz", "-o", str(filt_vcf), str(filt_vcf_plain)], log=log)
        index_vcf(filt_vcf, threads, log)

    ckpt.run(
        f"lea.maf_missing_filter.{prefix}",
        inputs={"biallelic_vcf": bial},
        outputs=[filt_vcf, Path(str(filt_vcf) + ".csi")],
        params={"threads": threads, "maf": maf, "max_missing": max_missing},
        action=action_vcftools_filter,
    )

    def action_plink_bed() -> None:
        run(["plink", "--vcf", str(filt_vcf), "--double-id", "--allow-extra-chr", "--make-bed", "--out", str(plink_prefix)], log=log)

    ckpt.run(
        f"lea.plink_bed.{prefix}",
        inputs={"filtered_vcf": filt_vcf},
        outputs=bed_outputs(plink_prefix),
        params={"double_id": True},
        action=action_plink_bed,
    )

    prune_in = Path(str(prune_prefix) + ".prune.in")

    def action_prune_list() -> None:
        run(["plink", "--bfile", str(plink_prefix), "--allow-extra-chr", "--indep-pairwise", window, step, r2, "--out", str(prune_prefix)], log=log)

    ckpt.run(
        f"lea.prune_list.{prefix}",
        inputs={"bed": Path(str(plink_prefix) + ".bed"), "bim": Path(str(plink_prefix) + ".bim"), "fam": Path(str(plink_prefix) + ".fam")},
        outputs=[prune_in],
        params={"window": window, "step": step, "r2": r2},
        action=action_prune_list,
    )

    def action_extract_pruned() -> None:
        run(["plink", "--bfile", str(plink_prefix), "--extract", str(prune_in), "--allow-extra-chr", "--make-bed", "--out", str(lea_bed)], log=log)
        run(["plink", "--bfile", str(lea_bed), "--allow-extra-chr", "--recode", "--out", str(lea_final)], log=log)
        run(["bcftools", "view", "--threads", str(threads), "-i", f"ID=@{prune_in}", str(filt_vcf), "-Oz", "-o", str(pruned_vcf)], log=log)
        index_vcf(pruned_vcf, threads, log)

    ckpt.run(
        f"lea.extract_pruned.{prefix}",
        inputs={"plink_bed": Path(str(plink_prefix) + ".bed"), "plink_bim": Path(str(plink_prefix) + ".bim"), "plink_fam": Path(str(plink_prefix) + ".fam"), "prune_in": prune_in, "filtered_vcf": filt_vcf},
        outputs=bed_outputs(lea_bed) + [Path(str(lea_final) + ".ped"), Path(str(lea_final) + ".map"), pruned_vcf, Path(str(pruned_vcf) + ".csi")],
        params={"extract_pruned": True},
        action=action_extract_pruned,
    )

    if cfg.get("_phase") == "prepare":
        return

    if as_bool(cfg.get("run_snmf"), default=True):
        rscript = scripts_dir / "run_lea.R"
        lea_out = outdir / f"{prefix}.LEA"
        geno = Path(str(lea_out) + ".geno")
        project = Path(str(lea_out) + ".snmfProject")
        ce = Path(str(lea_out) + ".cross_entropy.tsv")
        ce_summary = Path(str(lea_out) + ".cross_entropy_summary.tsv")
        summary_done = Path(str(lea_out) + ".summary.done")

        k_min = str(cfg.get("k_min", 2))
        k_max = str(cfg.get("k_max", 20))
        repetitions = str(cfg.get("repetitions", 5))
        seed = str(cfg.get("seed", 12345))

        def base_rargs(extra: Sequence[str]) -> List[str]:
            args = [
                "Rscript", str(rscript),
                "--ped", str(lea_final) + ".ped",
                "--fam", str(lea_bed) + ".fam",
                "--out", str(lea_out),
                "--k-min", k_min,
                "--k-max", k_max,
                "--repetitions", repetitions,
                "--seed", seed,
                "--threads", str(threads),
                "--plot-k-min", str(cfg.get("plot_k_min", cfg.get("k_min", 2))),
                "--plot-k-max", str(cfg.get("plot_k_max", cfg.get("k_max", 20))),
                "--order-by", str(cfg.get("order_by", "bestK")),
                "--width", str(cfg.get("width", 20)),
                "--height", str(cfg.get("height", 6)),
                "--multi-height", str(cfg.get("multi_height", 18)),
            ]
            maybe_add_path_arg(args, "--rename", rename)
            maybe_add_path_arg(args, "--clades", clades)
            maybe_add_path_arg(args, "--clade_colors", clade_colors)
            if color_sample_labels_by_clade:
                args.append("--color_sample_labels_by_clade")
            if clade_color_strip:
                args.append("--clade_color_strip")
            return args + list(extra)

        def action_snmf_model() -> None:
            run(base_rargs(["--model-only"]), log=log)

        model_name = f"lea.snmf_model.{prefix}"
        model_inputs = {"ped": Path(str(lea_final) + ".ped"), "fam": Path(str(lea_bed) + ".fam"), "rscript": rscript}
        model_outputs = [geno, project]
        model_params = {"k_min": k_min, "k_max": k_max, "repetitions": repetitions, "seed": seed, "threads": threads}
        if cfg.get("_phase") != "plot":
            ckpt.run(
                model_name,
                inputs=model_inputs,
                outputs=model_outputs,
                params=model_params,
                action=action_snmf_model,
            )
        if cfg.get("_phase") == "model":
            return

        best_k = cfg.get("best_k", cfg.get("K", cfg.get("k", "auto")))
        plot_k_min = str(cfg.get("plot_k_min", cfg.get("k_min", 2)))
        plot_k_max = str(cfg.get("plot_k_max", cfg.get("k_max", 20)))
        plot_outputs = [
            ce,
            ce_summary,
            Path(str(lea_out) + ".bestK.txt"),
            Path(str(lea_out) + ".best.Q.tsv"),
            Path(str(lea_out) + ".best.admixture.pdf"),
            Path(str(lea_out) + ".cross_entropy.pdf"),
            Path(str(lea_out) + f".K{plot_k_min}_to_K{plot_k_max}.admixture.pdf"),
            summary_done,
        ]
        if clades is not None:
            plot_outputs += [
                Path(str(lea_out) + ".clade_colors.tsv"),
                Path(str(lea_out) + ".sample_clades.tsv"),
            ]
        if best_k is not None and str(best_k).strip().lower() != "auto":
            plot_outputs += [
                Path(str(lea_out) + f".K{str(best_k)}.Q.tsv"),
                Path(str(lea_out) + f".K{str(best_k)}.admixture.pdf"),
            ]

        def action_snmf_plots() -> None:
            extra = ["--summary-only"]
            if best_k is not None:
                extra.extend(["--best-k", str(best_k)])
            run(base_rargs(extra), log=log)

        ckpt.run(
            f"lea.snmf_summarize_plot.{prefix}",
            inputs={"geno": geno, "project": project, "fam": Path(str(lea_bed) + ".fam"), "rename": rename, "clades": clades, "clade_colors": clade_colors, "rscript": rscript, "palette_helper": scripts_dir / "clade_palette_utils.R"},
            outputs=plot_outputs,
            params={"best_k": best_k, "plot_k_min": plot_k_min, "plot_k_max": plot_k_max, "order_by": str(cfg.get("order_by", "bestK")), "width": cfg.get("width", 20), "height": cfg.get("height", 6), "multi_height": cfg.get("multi_height", 18), "rename": str(rename) if rename else None, "clades": str(clades) if clades else None, "clade_colors": str(clade_colors) if clade_colors else None, "color_sample_labels_by_clade": color_sample_labels_by_clade, "clade_color_strip": clade_color_strip},
            action=action_snmf_plots,
        )

    if as_bool(cfg.get("cleanup_intermediates"), default=False):
        print("[lea] cleanup_intermediates=true, but checkpointed intermediates are retained for reuse")
