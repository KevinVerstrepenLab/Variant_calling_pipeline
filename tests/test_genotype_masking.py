"""Verify sample-specific optional masks with the actual bcftools plugin."""
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow'))
BCFTOOLS = os.environ.get('BCFTOOLS_BINARY') or shutil.which('bcftools')
HAS_SNAKEMAKE = importlib.util.find_spec('snakemake') is not None


@unittest.skipUnless(BCFTOOLS and HAS_SNAKEMAKE, 'bcftools and Snakemake required')
class MaskingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.vcf = Path(self.temp.name) / 'input.vcf'
        self.vcf.write_text('''##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tref\tmixed\talt\tlow
chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT:DP:AD\t0/0:10:10,0\t0/1:10:5,5\t1/1:10:0,10\t0/0:4:4,0
''')

    def mask(self, expression):
        prefix = Path(BCFTOOLS).resolve().parents[1]
        env = dict(os.environ, BCFTOOLS_PLUGINS=str(prefix / 'libexec/bcftools'))
        result = subprocess.run([BCFTOOLS, '+setGT', str(self.vcf), '-Ov', '--',
                                 '-t', 'q', '-n', '.', '-i', expression],
                                env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        record = next(line for line in result.stdout.splitlines() if not line.startswith('#'))
        return [sample.split(':')[0] for sample in record.split('\t')[9:]]

    def test_ad_and_low_depth_mask_only_affected_samples(self):
        from common import ad_expr
        self.assertEqual(self.mask(ad_expr({'ad_filter': {'min_dp': 5}})),
                         ['0/0', './.', '1/1', './.'])

    def test_heterozygote_mask_preserves_homozygotes(self):
        self.assertEqual(self.mask('GT="het"'), ['0/0', './.', '1/1', '0/0'])


if __name__ == '__main__':
    unittest.main()
