"""Translate the sectioned run config for existing rules and validate analysis settings."""
from copy import deepcopy
from pathlib import Path
import re


def normalize_config(raw):
    """Keep legacy flat keys working; explicit legacy keys override section defaults."""
    cfg = deepcopy(raw)
    sections = {
        "inputs": {
            "reference": "ref", "samples_tsv": "samples_tsv", "fastq_dir": "fastq_dir",
            "read_layout": "read_layout", "single_end_pattern": "single_end_pattern",
        },
        "project": {"outdir": "outdir"},
        "read_processing": {"call_bam": "call_bam"},
        "variant_calling": {"gatk": "gatk", "joint": "joint", "gvcf_reuse": "gvcf_reuse"},
        "filtering": {"hard": "filters", "genotypes": "ad_filter", "final": "final"},
    }
    for section, keys in sections.items():
        values = cfg.get(section, {})
        if not isinstance(values, dict):
            raise ValueError(f"{section} must be a mapping")
        for new, old in keys.items():
            if new in values:
                cfg.setdefault(old, values[new])
    for stage, settings in cfg.get("resources", {}).items():
        if not isinstance(settings, dict):
            raise ValueError(f"resources.{stage} must be a mapping")
        for key in ("threads", "mem_mb"):
            if key in settings:
                cfg.setdefault(f"{stage}_{key}", settings[key])
    for key, value in cfg.get("read_processing", {}).get("fastp", {}).items():
        cfg.setdefault(f"fastp_{key}", value)
    return cfg


def postprocessing_settings(config):
    """Resolve shared metadata once, using the workflow working directory for paths."""
    pp = deepcopy(config.get("postprocessing", {}))
    pp.setdefault("enabled", False)
    pp.setdefault("prefix", "cohort")
    # PCA and trees share variant preparation, but have separate user settings.
    if any(key in pp for key in ("variant_selection", "pca", "phylogenetic_tree")):
        pca = pp.get("pca", {})
        phylogeny = pp.get("phylogenetic_tree", {})
        combined = deepcopy(pp.get("variant_selection", {}))
        combined.update({key: value for key, value in phylogeny.items() if key not in {"enabled", "sample_list"}})
        combined.update({
            "enabled": pca.get("enabled", True) or phylogeny.get("enabled", True),
            "pca": pca.get("enabled", True),
            "tree": phylogeny.get("enabled", True),
            "pca_components": pca.get("components", 20),
            "pca_label_points": pca.get("label_points", False),
            "pca_label_size": pca.get("label_size", 2.0),
            "pca_sample_list": pca.get("sample_list"),
            "tree_sample_list": phylogeny.get("sample_list"),
        })
        # Preserve explicit overrides in earlier combined configurations.
        combined.update(pp.get("pca_tree", {}))
        pp["pca_tree"] = combined
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(pp["prefix"])):
        raise ValueError("postprocessing.prefix must contain only letters, digits, _, . or -")
    metadata = pp.get("metadata", {}) or {}
    for stage in ("pca_tree", "lea"):
        settings = pp.setdefault(stage, {})
        settings.setdefault("enabled", True)
        settings["force"] = pp.get("force", False)
        for key in ("rename", "clades", "clade_colors"):
            if settings.get(key) is None:
                settings[key] = metadata.get(key)
        for key in ("keep", "rename", "clades", "clade_colors", "pca_sample_list", "tree_sample_list"):
            if settings.get(key):
                path = Path(settings[key]).expanduser().resolve()
                if pp["enabled"] and settings["enabled"] and not path.is_file():
                    raise ValueError(f"postprocessing.{stage}.{key} not found: {path}")
                if pp["enabled"] and settings["enabled"] and path.is_file() and path.stat().st_size == 0:
                    raise ValueError(f"postprocessing.{stage}.{key} is empty: {path}")
                settings[key] = str(path)
    tree = pp["pca_tree"]
    if pp["enabled"] and tree["enabled"]:
        if not (tree.get("pca", True) or tree.get("tree", True)):
            raise ValueError("Enable pca or tree, or disable postprocessing.pca_tree")
        if tree.get("method", "nj") not in {"nj", "upgma"}:
            raise ValueError("postprocessing.phylogenetic_tree.method must be nj or upgma")
        if not re.fullmatch(r"[A-Za-z_]+", tree.get("layout", "fan")):
            raise ValueError("Invalid postprocessing.phylogenetic_tree.layout")
        if tree.get("pca", True) and int(tree.get("pca_components", 20)) < 2:
            raise ValueError("postprocessing.pca.components must be at least 2")
    lea = pp["lea"]
    if pp["enabled"] and lea["enabled"] and lea.get("run_snmf", True):
        lo, hi = int(lea.get("k_min", 2)), int(lea.get("k_max", 20))
        if lo < 1 or hi < lo or int(lea.get("repetitions", 5)) < 1:
            raise ValueError("LEA requires 1 <= k_min <= k_max and repetitions >= 1")
        best = lea.get("best_k", "auto")
        if str(best).lower() != "auto" and not lo <= int(best) <= hi:
            raise ValueError("postprocessing.lea.best_k must be auto or inside the tested K range")
        # Use the same clipped range in Python's output contract and the R plotter.
        lea["plot_k_min"] = max(lo, int(lea.get("plot_k_min", lo)))
        lea["plot_k_max"] = min(hi, int(lea.get("plot_k_max", hi)))
        if lea["plot_k_min"] > lea["plot_k_max"]:
            raise ValueError("LEA plot K range must overlap the tested K range")
    return pp


def postprocessing_phase_settings(settings, stage, phase):
    """Separate preparation/model parameters from options that only affect plots."""
    cfg = settings[stage]
    shared = {"enabled", "force", "maf", "ld_prune"}
    if stage == "pca_tree":
        shared |= {"geno", "snps_only", "biallelic", "keep", "pca", "tree", "pca_components"}
    else:
        shared |= {"max_missing", "run_snmf"}
        if phase == "model":
            shared |= {"k_min", "k_max", "repetitions", "seed"}
    if phase == "plot":
        return deepcopy(cfg)
    return {key: deepcopy(value) for key, value in cfg.items() if key in shared}


def postprocessing_outputs(settings, outdir, stage, phase=None):
    """Concrete products for each native Snakemake analysis rule."""
    cfg = settings[stage]
    if not settings["enabled"] or not cfg["enabled"]:
        return []
    prefix = f"{outdir}/postprocessing/{stage}/{settings['prefix']}"
    prepare, model, plot = [], [], []
    if stage == "pca_tree":
        if cfg.get("pca", True):
            prepare += [prefix + ".pca" + ext for ext in (".eigenvec", ".eigenval")]
            plot += [prefix + ".tree" + ext for ext in (".pca_scores.tsv", ".pca.pdf", ".pca.png")]
        if cfg.get("tree", True):
            prepare += [prefix + ".dist" + ext for ext in (".mdist", ".mdist.id")]
            plot += [prefix + ".tree" + ext for ext in (".nwk", ".tree_info.tsv")]
            plot += [prefix + f".tree-{cfg.get('layout', 'fan')}" + ext for ext in (".pdf", ".png")]
        if cfg.get("clades"):
            plot.append(prefix + ".tree.clade_colors.tsv")
    else:
        prepare = [prefix + ".LEA.final" + ext for ext in (".ped", ".map")]
        prepare += [prefix + ".LEA" + ext for ext in (".bed", ".bim", ".fam")]
        prepare += [prefix + ".pruned.vcf.gz", prefix + ".pruned.vcf.gz.csi"]
        if cfg.get("run_snmf", True):
            model = [prefix + ".LEA" + ext for ext in (".geno", ".snmfProject")]
            for k in range(int(cfg.get("k_min", 2)), int(cfg.get("k_max", 20)) + 1):
                for rep in range(1, int(cfg.get("repetitions", 5)) + 1):
                    model += [prefix + f".LEA.snmf/K{k}/run{rep}/{settings['prefix']}.LEA_r{rep}.{k}" + ext
                              for ext in (".Q", ".G")]
            plot = [prefix + ".LEA" + ext for ext in (
                ".cross_entropy.tsv", ".cross_entropy_summary.tsv", ".bestK.txt",
                ".best.Q.tsv", ".best.admixture.pdf", ".cross_entropy.pdf", ".summary.done",
                f".K{cfg['plot_k_min']}_to_K{cfg['plot_k_max']}.admixture.pdf",
            )]
            if cfg.get("clades"):
                plot += [prefix + ".LEA" + ext for ext in (".clade_colors.tsv", ".sample_clades.tsv")]
    products = {"prepare": prepare, "model": model, "plot": plot}
    return products[phase] if phase else prepare + model + plot


def coverage_settings(config):
    """Defaults and early validation for optional local coverage screening."""
    settings = {
        'enabled': False, 'bam_dir': None, 'window_bp': 10000, 'mapq': 30,
        'exclude_flags': 3844, 'min_contig_bp': 10000, 'exclude_contigs': [],
        'baseline_contigs': [], 'max_n_fraction': 0.10, 'gc_correct': True,
        'gc_bins': 20, 'min_gc_windows': 200, 'smooth_windows': 5,
        'log2_threshold': 0.40, 'min_region_windows': 3, 'pseudocount': 0.5,
        'min_baseline_depth': 5, 'width': 18, 'height': 8, 'dpi': 200,
        'plot_min_contig_bp': None, 'chromosomes_per_row': 16,
        'outlier_quantiles': [0.01, 0.99],
    }
    settings.update(config.get('coverage_plots', {}) or {})
    if settings['plot_min_contig_bp'] is None:
        settings['plot_min_contig_bp'] = 50000
    if settings['chromosomes_per_row'] is None:
        settings['chromosomes_per_row'] = 16
    value = settings['plot_min_contig_bp']
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError('coverage_plots.plot_min_contig_bp must be null or a nonnegative integer')
    value = settings['chromosomes_per_row']
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError('coverage_plots.chromosomes_per_row must be a positive integer')
    quantiles = settings['outlier_quantiles']
    if quantiles is not None and (not isinstance(quantiles, list) or len(quantiles) != 2
            or any(isinstance(q, bool) or not isinstance(q, (int, float)) for q in quantiles)
            or not 0 <= quantiles[0] < quantiles[1] <= 1):
        raise ValueError('coverage_plots.outlier_quantiles must be null or [lower, upper] between 0 and 1')
    for key in ('enabled', 'gc_correct'):
        if not isinstance(settings[key], bool):
            raise ValueError(f'coverage_plots.{key} must be true or false')
    for key in ('window_bp', 'min_contig_bp', 'gc_bins', 'min_gc_windows',
                'smooth_windows', 'min_region_windows', 'dpi'):
        value = settings[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f'coverage_plots.{key} must be a positive integer')
    if settings['smooth_windows'] % 2 != 1:
        raise ValueError('coverage_plots.smooth_windows must be odd')
    if settings['gc_bins'] < 2 or settings['min_gc_windows'] < 20:
        raise ValueError('coverage_plots requires gc_bins >= 2 and min_gc_windows >= 20')
    for key in ('mapq', 'exclude_flags'):
        if not isinstance(settings[key], int) or settings[key] < 0:
            raise ValueError(f'coverage_plots.{key} must be a nonnegative integer')
    if not 0 <= settings['max_n_fraction'] <= 1:
        raise ValueError('coverage_plots.max_n_fraction must be between 0 and 1')
    for key in ('log2_threshold', 'pseudocount', 'width', 'height'):
        if settings[key] <= 0:
            raise ValueError(f'coverage_plots.{key} must be positive')
    if settings['min_baseline_depth'] < 0:
        raise ValueError('coverage_plots.min_baseline_depth must be nonnegative')
    for key in ('exclude_contigs', 'baseline_contigs'):
        if not isinstance(settings[key], list) or any(not isinstance(value, str) for value in settings[key]):
            raise ValueError(f'coverage_plots.{key} must be a list of exact contig names')
    if settings['bam_dir']:
        settings['bam_dir'] = str(Path(settings['bam_dir']).expanduser().resolve())
        if settings['enabled'] and not Path(settings['bam_dir']).is_dir():
            raise ValueError('coverage_plots.bam_dir does not exist')
    return settings
