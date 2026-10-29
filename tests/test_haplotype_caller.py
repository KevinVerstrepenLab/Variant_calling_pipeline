"""Run the real HaplotypeCaller rule with a fake GATK to test publication failures.

Run with the project's Snakemake environment, for example:
  python -m pytest -q tests/test_haplotype_caller.py

Set HC_RULE_FILE to validate a candidate copy of the production rule before it
is installed. Fixtures and all Snakemake state live under /tmp.
"""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap
import unittest

ROOT = Path(os.environ.get('HC_TEST_ROOT', Path(__file__).resolve().parents[1]))
RULE = Path(os.environ.get('HC_RULE_FILE', ROOT / 'workflow/rules/haplotype_caller.smk'))
SNAKEMAKE = os.environ.get('SNAKEMAKE_BINARY') or shutil.which('snakemake')


@unittest.skipUnless(SNAKEMAKE, 'Snakemake executable required')
class HaplotypeCallerPublicationTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix='hc-publication-', dir='/tmp')
        self.addCleanup(self.scratch.cleanup)
        self.work = Path(self.scratch.name) / "run with 'quoted paths'"
        self.work.mkdir()
        self.gvcf_dir = self.work / 'gvcf'
        self.gvcf_dir.mkdir()
        self.final = self.gvcf_dir / 'sample.g.vcf.gz'
        self.trace = self.work / 'gatk-arguments.json'
        self.tmpdir = self.work / 'java temporary files'
        self.tmpdir.mkdir()
        self.fake_bin = self.work / 'bin'
        self.fake_bin.mkdir()
        fixture_dir = self.work / 'inputs'
        fixture_dir.mkdir()
        for filename in ('sample.call.bam', 'sample.call.bam.bai', 'reference.fa',
                         'reference.fa.fai', 'reference.dict'):
            (fixture_dir / filename).write_text('test fixture\n')

        (self.work / 'Snakefile').write_text(textwrap.dedent(f'''\
            from pathlib import Path
            import re
            REF = {str(fixture_dir / 'reference.fa')!r}
            GVCF = {str(self.gvcf_dir)!r}
            LOGS = {str(self.work / 'logs')!r}
            ENVS = Path({str(ROOT / 'workflow/envs')!r})
            TMP_DIR = {str(self.tmpdir)!r}
            CALL_SAMPLES = ['sample']
            config = {{'hc_threads': 1, 'hc_mem_mb': 16000}}
            rule make_call_bam:
                output:
                    bam={str(fixture_dir / '{sample}.call.bam')!r},
                    bai={str(fixture_dir / '{sample}.call.bam.bai')!r}
            rule ref_fai:
                output: {str(fixture_dir / 'reference.fa.fai')!r}
            rule ref_dict:
                output: {str(fixture_dir / 'reference.dict')!r}
            include: {str(RULE)!r}
        '''))
        fake_gatk = self.fake_bin / 'gatk'
        fake_gatk.write_text('#!/usr/bin/env python3\n' + textwrap.dedent('''\
            import gzip
            import json
            import os
            from pathlib import Path
            import sys
            args = sys.argv[1:]
            Path(os.environ['HC_TRACE']).write_text(json.dumps(args))
            output = Path(args[args.index('-O') + 1])
            mode = os.environ['HC_FAKE_MODE']
            if mode == 'partial-failure':
                output.write_bytes(b'partial, invalid compressed GVCF')
                sys.exit(3)
            if mode == 'empty-output':
                output.write_bytes(b'')
            else:
                with gzip.open(output, 'wt') as handle:
                    handle.write('##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tsample\\n')
            if mode not in ('missing-index', 'empty-index'):
                Path(str(output) + '.tbi').write_bytes(b'complete index')
            elif mode == 'empty-index':
                Path(str(output) + '.tbi').write_bytes(b'')
        '''))
        fake_gatk.chmod(0o755)

    def run_rule(self, mode):
        env = dict(os.environ,
                   PATH=str(self.fake_bin) + os.pathsep + os.environ['PATH'],
                   HC_TRACE=str(self.trace), HC_FAKE_MODE=mode,
                   XDG_CACHE_HOME=str(self.work / 'cache'))
        completed = subprocess.run([
            SNAKEMAKE, str(self.final), '--snakefile', str(self.work / 'Snakefile'),
            '--cores', '1', '--printshellcmds', '--nolock',
            '--allowed-rules', 'haplotype_caller_gvcf',
        ], cwd=self.work, env=env, capture_output=True, text=True, timeout=60)
        self.last_log = completed.stdout + completed.stderr
        return completed

    def assert_staging_cleaned(self):
        self.assertEqual(list(self.gvcf_dir.glob('.sample.hc.*')), [], self.last_log)

    def assert_gatk_used_staging(self):
        arguments = json.loads(self.trace.read_text())
        staged_output = Path(arguments[arguments.index('-O') + 1])
        self.assertNotEqual(staged_output, self.final)
        self.assertEqual(staged_output.parent.parent, self.gvcf_dir)
        self.assertTrue(staged_output.parent.name.startswith('.sample.hc.'))
        return arguments

    def assert_failed_output_unpublished(self, mode):
        completed = self.run_rule(mode)
        self.assertNotEqual(completed.returncode, 0, self.last_log)
        self.assertTrue(self.trace.is_file(), 'fake GATK was not invoked:\n' + self.last_log)
        self.assert_gatk_used_staging()
        self.assertEqual(list(self.gvcf_dir.glob('*.g.vcf.gz')), [], self.last_log)
        self.assertFalse(Path(str(self.final) + '.tbi').exists(), self.last_log)
        self.assert_staging_cleaned()

    def test_partial_output_is_not_published_on_gatk_failure(self):
        self.assert_failed_output_unpublished('partial-failure')

    def test_success_publishes_complete_gvcf_and_index(self):
        completed = self.run_rule('success')
        self.assertEqual(completed.returncode, 0, self.last_log)
        self.assertGreater(self.final.stat().st_size, 0)
        self.assertGreater(Path(str(self.final) + '.tbi').stat().st_size, 0)
        with gzip.open(self.final, 'rt') as handle:
            self.assertIn('##fileformat=VCFv4.2', handle.read())
        arguments = self.assert_gatk_used_staging()
        self.assertIn('--tmp-dir', arguments)
        self.assertEqual(arguments[arguments.index('--tmp-dir') + 1], str(self.tmpdir))
        self.assert_staging_cleaned()

    def test_missing_index_is_not_published(self):
        self.assert_failed_output_unpublished('missing-index')

    def test_empty_index_is_not_published(self):
        self.assert_failed_output_unpublished('empty-index')

    def test_empty_gvcf_is_not_published(self):
        self.assert_failed_output_unpublished('empty-output')


if __name__ == '__main__':
    unittest.main()
