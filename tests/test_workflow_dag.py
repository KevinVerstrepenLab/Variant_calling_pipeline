"""Exercise mixed layouts, GVCF reuse and cohort selection with Snakemake."""
import gzip
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = os.environ.get('SNAKEMAKE_BINARY') or shutil.which('snakemake')


@unittest.skipUnless(SNAKEMAKE, 'Snakemake executable required')
class DagTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out = Path(self.temp.name) / 'results'
        (self.out / 'gvcf').mkdir(parents=True)

    def existing(self, sample):
        path = self.out / 'gvcf' / f'{sample}.g.vcf.gz'
        with gzip.open(path, 'wt') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n')
        return path

    def dry_run(self, include_existing=False):
        env = dict(os.environ, XDG_CACHE_HOME=str(Path(self.temp.name) / 'cache'))
        config = Path(self.temp.name) / 'override.json'
        config.write_text(json.dumps({
            'outdir': str(self.out),
            'gvcf_reuse': {'enabled': True, 'include_existing_samples': include_existing},
        }))
        completed = subprocess.run([
            SNAKEMAKE, '--snakefile', 'workflow/Snakefile', '--configfile',
            'tests/minimal/config.yaml', str(config), '--dry-run', '--cores', '2',
        ], cwd=ROOT, env=env, capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        return completed.stdout + completed.stderr

    def test_mixed_layouts_have_both_trimming_rules(self):
        output = self.dry_run()
        self.assertIn('fastp_paired', output)
        self.assertIn('fastp_single', output)

    def test_reused_sample_has_no_upstream_jobs(self):
        path = self.existing('paired')
        before = path.read_bytes()
        output = self.dry_run()
        self.assertNotIn('fastp_paired', output)
        self.assertIn('fastp_single', output)
        self.assertEqual(path.read_bytes(), before)

    def test_extra_samples_require_opt_in(self):
        self.existing('unlisted')
        self.assertNotIn('unlisted.g.vcf.gz', self.dry_run())
        self.assertIn('unlisted.g.vcf.gz', self.dry_run(include_existing=True))

    def test_plain_legacy_sample_is_reused_with_tribble_index(self):
        path = self.existing('paired')
        plain = path.with_name('paired.snps.indels.g.vcf')
        with gzip.open(path, 'rb') as source:
            plain.write_bytes(source.read())
        path.unlink()
        output = self.dry_run()
        self.assertNotIn('fastp_paired', output)
        self.assertIn('paired.snps.indels.g.vcf.idx', output)
        self.assertIn('index_legacy_gvcf', output)


if __name__ == '__main__':
    unittest.main()
