"""Regression checks for the unified config, analysis DAG and model reuse."""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow"))
from configuration import normalize_config, postprocessing_outputs, postprocessing_settings, postprocessing_phase_settings

spec = importlib.util.spec_from_file_location("postprocess_vcf", ROOT / "workflow/scripts/postprocess_vcf.py")
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
SNAKEMAKE = os.environ.get("SNAKEMAKE_BINARY") or shutil.which("snakemake")


class ConfigTests(unittest.TestCase):
    def test_sample_lists_are_resolved_for_plotting_only(self):
        with tempfile.TemporaryDirectory() as temp:
            pca = Path(temp) / "pca.txt"
            tree = Path(temp) / "tree.txt"
            pca.write_text("sample1\n")
            tree.write_text("sample2\nsample3\n")
            settings = postprocessing_settings({"postprocessing": {
                "enabled": True, "pca": {"sample_list": str(pca)},
                "phylogenetic_tree": {"sample_list": str(tree)},
            }})
            plot = postprocessing_phase_settings(settings, "pca_tree", "plot")
            prepare = postprocessing_phase_settings(settings, "pca_tree", "prepare")
            self.assertEqual(plot["pca_sample_list"], str(pca.resolve()))
            self.assertEqual(plot["tree_sample_list"], str(tree.resolve()))
            self.assertNotIn("pca_sample_list", prepare)
            self.assertNotIn("tree_sample_list", prepare)

    def test_separate_pca_and_tree_sections(self):
        settings = postprocessing_settings({"postprocessing": {
            "enabled": True,
            "variant_selection": {"maf": 0.03, "geno": 0.2},
            "pca": {"enabled": True, "components": 5, "label_points": True},
            "phylogenetic_tree": {"enabled": False, "method": "upgma", "label_tips": False},
            "lea": {"enabled": False},
        }})
        combined = settings["pca_tree"]
        self.assertTrue(combined["pca"])
        self.assertFalse(combined["tree"])
        self.assertEqual(combined["pca_components"], 5)
        self.assertEqual(combined["maf"], 0.03)
        self.assertTrue(combined["pca_label_points"])
        self.assertFalse(combined["label_tips"])
        self.assertFalse(any(path.endswith(".nwk") for path in postprocessing_outputs(settings, "results", "pca_tree")))

    def test_both_pca_and_tree_disabled(self):
        settings = postprocessing_settings({"postprocessing": {
            "enabled": True, "pca": {"enabled": False}, "phylogenetic_tree": {"enabled": False},
        }})
        self.assertEqual(postprocessing_outputs(settings, "results", "pca_tree"), [])

    def test_nested_config_and_legacy_overrides(self):
        cfg = normalize_config({
            "inputs": {"reference": "nested.fa"}, "ref": "legacy.fa",
            "project": {"outdir": "new-results"},
            "resources": {"fastp": {"threads": 3, "mem_mb": 5000}},
            "read_processing": {"fastp": {"min_len": 40}},
            "filtering": {"final": {"max_missing": 0.9}},
        })
        self.assertEqual(cfg["ref"], "legacy.fa")
        self.assertEqual(cfg["outdir"], "new-results")
        self.assertEqual(cfg["fastp_threads"], 3)
        self.assertEqual(cfg["fastp_mem_mb"], 5000)
        self.assertEqual(cfg["fastp_min_len"], 40)
        self.assertEqual(cfg["final"]["max_missing"], 0.9)

    def test_pca_only_output_contract(self):
        settings = postprocessing_settings({"postprocessing": {
            "enabled": True, "pca_tree": {"tree": False}, "lea": {"enabled": False},
        }})
        outputs = postprocessing_outputs(settings, "results", "pca_tree")
        self.assertTrue(any(path.endswith(".pca.pdf") for path in outputs))
        self.assertFalse(any(path.endswith(".nwk") for path in outputs))
        self.assertEqual(postprocessing_outputs(settings, "results", "lea"), [])

    def test_clipped_plot_range_and_stable_best_products(self):
        settings = postprocessing_settings({"postprocessing": {
            "enabled": True, "lea": {"k_min": 3, "k_max": 5, "plot_k_min": 2, "plot_k_max": 12},
        }})
        outputs = postprocessing_outputs(settings, "results", "lea")
        self.assertTrue(any(path.endswith(".K3_to_K5.admixture.pdf") for path in outputs))
        self.assertTrue(any(path.endswith(".best.Q.tsv") for path in outputs))

    def test_invalid_best_k_and_missing_metadata_fail_early(self):
        for settings in ({"lea": {"k_max": 3, "best_k": 8}},
                         {"metadata": {"clades": "/does/not/exist.tsv"}}):
            with self.assertRaises(ValueError):
                postprocessing_settings({"postprocessing": {"enabled": True, **settings}})

    def test_plot_parameters_are_excluded_from_model_and_preparation_jobs(self):
        settings = postprocessing_settings({"postprocessing": {
            "enabled": True, "lea": {"seed": 19, "best_k": 4, "order_by": "clade"},
        }})
        preparation = postprocessing_phase_settings(settings, "lea", "prepare")
        model = postprocessing_phase_settings(settings, "lea", "model")
        self.assertNotIn("seed", preparation)
        self.assertEqual(model["seed"], 19)
        self.assertNotIn("best_k", model)
        self.assertNotIn("order_by", model)
        self.assertNotIn("clades", model)
        self.assertEqual(postprocessing_phase_settings(settings, "lea", "plot")["best_k"], 4)


class CheckpointTests(unittest.TestCase):
    def test_model_changes_rerun_but_plot_changes_reuse_model(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ped = root / "input.ped"
            ped.write_text("original genotypes")
            model = root / "model"
            checkpoint = engine.Checkpointer(root / "checkpoints")
            executions = []

            def model_action():
                executions.append(True)
                model.write_text("model output")

            def execute(seed):
                with contextlib.redirect_stdout(io.StringIO()):
                    checkpoint.run("model", {"ped": ped}, [model], {"seed": seed}, model_action)

            execute(123)
            execute(123)
            self.assertEqual(len(executions), 1)
            execute(456)
            self.assertEqual(len(executions), 2)
            ped.write_text("changed genotypes")
            execute(456)
            self.assertEqual(len(executions), 3)


@unittest.skipUnless(SNAKEMAKE, "Snakemake executable required")
class AnalysisDagTests(unittest.TestCase):
    def test_both_analyses_use_final_snps_and_shared_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            clades = root / "clades.tsv"
            clades.write_text("sample\tclade\npaired\tA\n")
            cfg = root / "config.json"
            cfg.write_text(json.dumps({
                "inputs": {"reference": "tests/minimal/reference.fa", "samples_tsv": "tests/minimal/samples.tsv"},
                "project": {"outdir": str(root / "results")},
                "postprocessing": {"enabled": True, "metadata": {"clades": str(clades)}},
            }))
            result = subprocess.run([
                SNAKEMAKE, "--snakefile", "workflow/Snakefile", "--configfile", str(cfg),
                "--dry-run", "--cores", "2", "--printshellcmds",
            ], cwd=ROOT, env=dict(os.environ, XDG_CACHE_HOME=str(root / "cache")), capture_output=True, text=True)
            output = result.stdout + result.stderr
            self.assertEqual(result.returncode, 0, output)
            self.assertIn("postprocess_pca_tree", output)
            self.assertIn("postprocess_lea", output)
            for rule in ("postprocess_pca_tree_prepare", "postprocess_pca_tree_plot",
                         "postprocess_lea_prepare", "postprocess_lea_model", "postprocess_lea_plot"):
                self.assertIn(rule, output)
            self.assertIn("filtered_snps.final.vcf.gz", output)
            self.assertIn(str(clades), output)
            self.assertNotIn("postprocess_config.yaml", output)


if __name__ == "__main__":
    unittest.main()
