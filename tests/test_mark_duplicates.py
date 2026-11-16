"""Exercise the production MarkDuplicates rule with simulated tool failures."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
RULE = Path(os.environ.get('MARKDUP_RULE_FILE', ROOT / 'workflow/rules/mark_duplicates.smk'))
SNAKEMAKE = os.environ.get('SNAKEMAKE_BINARY') or shutil.which('snakemake')


@unittest.skipUnless(SNAKEMAKE, 'Snakemake executable required')
class MarkDuplicatesPublicationTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix='markdup-publication-')
        self.addCleanup(self.scratch.cleanup)
        self.work = Path(self.scratch.name) / "run with 'quoted paths'"
        self.work.mkdir()
        self.bam_dir = self.work / 'bam'
        self.bam_dir.mkdir()
        self.final = self.bam_dir / 'sample.dedup.bam'
        self.metrics = self.bam_dir / 'sample.dedup.metrics.txt'
        self.trace = self.work / 'gatk-arguments.json'
        self.samtools_trace = self.work / 'samtools-arguments.jsonl'
        self.tmpdir = self.work / 'java temporary files'
        self.tmpdir.mkdir()
        self.fake_bin = self.work / 'bin'
        self.fake_bin.mkdir()
        for filename in ('sample.sorted.bam', 'sample.sorted.bam.bai'):
            (self.bam_dir / filename).write_text('input fixture\n')

        (self.work / 'Snakefile').write_text(textwrap.dedent(f'''\
            from pathlib import Path
            BAM = {str(self.bam_dir)!r}
            LOGS = {str(self.work / 'logs')!r}
            ENVS = Path({str(ROOT / 'workflow/envs')!r})
            TMP_DIR = {str(self.tmpdir)!r}
            config = {{'markdup_threads': 1, 'markdup_mem_mb': 16000}}
            include: {str(RULE)!r}
        '''))
        self.make_executable('gatk', '''\
            import json
            import os
            from pathlib import Path
            import sys
            args = sys.argv[1:]
            Path(os.environ['MARKDUP_TRACE']).write_text(json.dumps(args))
            output = Path(args[args.index('-O') + 1])
            metrics = Path(args[args.index('-M') + 1])
            mode = os.environ['MARKDUP_FAKE_MODE']
            parts = Path(str(output) + '.parts')
            if parts.exists():
                sys.exit('Spark output parts already exists')
            parts.mkdir()
            (parts / 'part-00000').write_bytes(b'Spark output partition')
            output.write_bytes(b'partial BAM' if mode == 'gatk-failure' else b'complete BAM')
            if mode != 'missing-metrics':
                metrics.write_text('duplicate metrics\\n')
            if mode != 'no-sbi':
                Path(str(output) + '.sbi').write_bytes(b'Spark index')
            if mode == 'gatk-failure':
                sys.exit(3)
        ''')
        self.make_executable('samtools', '''\
            import json
            import os
            from pathlib import Path
            import sys
            args = sys.argv[1:]
            with Path(os.environ['MARKDUP_SAMTOOLS_TRACE']).open('a') as trace:
                trace.write(json.dumps(args) + '\\n')
            mode = os.environ['MARKDUP_FAKE_MODE']
            if args[0] == 'quickcheck':
                sys.exit(1 if mode == 'quickcheck-failure' else 0)
            if args[0] != 'index':
                sys.exit('unexpected samtools command')
            index = Path(args[-1])
            index.write_bytes(b'partial index' if mode == 'index-failure' else b'complete index')
            if mode == 'index-failure':
                sys.exit(2)
        ''')

    def make_executable(self, name, source):
        path = self.fake_bin / name
        path.write_text('#!/usr/bin/env python3\n' + textwrap.dedent(source))
        path.chmod(0o755)

    def run_rule(self, mode='success'):
        env = dict(os.environ,
                   PATH=str(self.fake_bin) + os.pathsep + os.environ['PATH'],
                   MARKDUP_TRACE=str(self.trace),
                   MARKDUP_SAMTOOLS_TRACE=str(self.samtools_trace),
                   MARKDUP_FAKE_MODE=mode,
                   XDG_CACHE_HOME=str(self.work / 'cache'))
        completed = subprocess.run([
            SNAKEMAKE, str(self.final), '--snakefile', str(self.work / 'Snakefile'),
            '--cores', '1', '--printshellcmds', '--nolock',
            '--allowed-rules', 'mark_duplicates',
        ], cwd=self.work, env=env, capture_output=True, text=True, timeout=60)
        self.last_log = completed.stdout + completed.stderr
        return completed

    def assert_staging_cleaned(self):
        self.assertEqual(list(self.bam_dir.glob('.sample.markdup.*')), [], self.last_log)

    def assert_gatk_used_staging(self):
        self.assertTrue(self.trace.is_file(), 'GATK was not invoked:\n' + self.last_log)
        arguments = json.loads(self.trace.read_text())
        staged_bam = Path(arguments[arguments.index('-O') + 1])
        staged_metrics = Path(arguments[arguments.index('-M') + 1])
        self.assertNotEqual(staged_bam, self.final)
        self.assertEqual(staged_bam.parent.parent, self.bam_dir)
        self.assertTrue(staged_bam.parent.name.startswith('.sample.markdup.'))
        self.assertEqual(staged_metrics.parent, staged_bam.parent)
        return arguments, staged_bam

    def assert_failed_outputs_unpublished(self, mode):
        completed = self.run_rule(mode)
        self.assertNotEqual(completed.returncode, 0, self.last_log)
        self.assert_gatk_used_staging()
        for output in (self.final, self.metrics, Path(str(self.final) + '.bai'),
                       Path(str(self.final) + '.sbi')):
            self.assertFalse(output.exists(), self.last_log)
        self.assert_staging_cleaned()

    def assert_successful_publication(self, mode='success'):
        completed = self.run_rule(mode)
        self.assertEqual(completed.returncode, 0, self.last_log)
        self.assertEqual(self.final.read_bytes(), b'complete BAM')
        self.assertEqual(Path(str(self.final) + '.bai').read_bytes(), b'complete index')
        self.assertEqual(self.metrics.read_text(), 'duplicate metrics\n')
        arguments, staged_bam = self.assert_gatk_used_staging()
        self.assertEqual(arguments[arguments.index('--tmp-dir') + 1], str(self.tmpdir))
        samtools_calls = [json.loads(line) for line in self.samtools_trace.read_text().splitlines()]
        self.assertEqual([call[0] for call in samtools_calls], ['quickcheck', 'index'])
        self.assertEqual(samtools_calls[0][-1], str(staged_bam))
        self.assertEqual(samtools_calls[1][-2:], [str(staged_bam), str(staged_bam) + '.bai'])
        self.assert_staging_cleaned()

    def test_success_publishes_bam_metrics_and_indexes(self):
        self.assert_successful_publication()
        self.assertEqual(Path(str(self.final) + '.sbi').read_bytes(), b'Spark index')

    def test_success_does_not_require_optional_spark_index(self):
        self.assert_successful_publication('no-sbi')
        self.assertFalse(Path(str(self.final) + '.sbi').exists())

    def test_existing_parts_do_not_block_retry_or_get_deleted(self):
        stale_parts = Path(str(self.final) + '.parts')
        stale_parts.mkdir()
        sentinel = stale_parts / 'previous-attempt'
        sentinel.write_text('preserved evidence\n')
        self.assert_successful_publication()
        self.assertEqual(sentinel.read_text(), 'preserved evidence\n')

    def test_gatk_failure_does_not_publish_partial_outputs(self):
        self.assert_failed_outputs_unpublished('gatk-failure')

    def test_index_failure_does_not_publish_partial_outputs(self):
        self.assert_failed_outputs_unpublished('index-failure')

    def test_quickcheck_failure_does_not_publish_partial_outputs(self):
        self.assert_failed_outputs_unpublished('quickcheck-failure')
        calls = [json.loads(line) for line in self.samtools_trace.read_text().splitlines()]
        self.assertEqual([call[0] for call in calls], ['quickcheck'])

    def test_missing_metrics_does_not_publish_partial_outputs(self):
        self.assert_failed_outputs_unpublished('missing-metrics')


if __name__ == '__main__':
    unittest.main()
