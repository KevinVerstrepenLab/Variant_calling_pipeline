"""Execute the real final-filter rules for a nonempty SNP VCF and empty INDEL VCF."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = os.environ.get('SNAKEMAKE_BINARY') or shutil.which('snakemake')
VCF_BIN = os.environ.get('VCFOPS_TEST_BIN')


@unittest.skipUnless(SNAKEMAKE and VCF_BIN, 'Set VCFOPS_TEST_BIN for real final-filter testing')
class FinalFilteringTests(unittest.TestCase):
    def test_empty_and_large_nonempty_vcfs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            reference = root/'ref.fa'
            reference.write_text('>chr1\n'+'A'*10000+'\n')
            out = root/'results'
            (out/'snps').mkdir(parents=True)
            (out/'indels').mkdir()
            header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n')
            env = dict(os.environ, PATH=str(Path(VCF_BIN))+os.pathsep+os.environ['PATH'],
                       XDG_CACHE_HOME=str(root/'cache'))
            for kind in ('snps', 'indels'):
                plain = root/f'{kind}.vcf'
                plain.write_text(header+(''.join(f'chr1\t{i}\t.\tA\tC\t60\tPASS\t.\tGT\t0/0\t0/1\n'
                                               for i in range(1, 5001)) if kind == 'snps' else ''))
                compressed = out/kind/f'filtered_{kind}.merged.vcf.gz'
                subprocess.run([str(Path(VCF_BIN)/'bcftools'), 'view', '-Oz', '-o', str(compressed), str(plain)],
                               env=env, check=True, capture_output=True)
                subprocess.run([str(Path(VCF_BIN)/'tabix'), '-p', 'vcf', str(compressed)], env=env, check=True)
            config = root/'config.json'
            config.write_text(json.dumps({'ref': str(reference), 'samples_tsv': 'tests/minimal/samples.tsv',
                                          'outdir': str(out)}))
            targets = [str(out/'final'/f'filtered_{kind}.final.vcf.gz') for kind in ('snps', 'indels')]
            result = subprocess.run([SNAKEMAKE, *targets, '--snakefile', 'workflow/Snakefile',
                '--configfile', str(config), '--cores', '2', '--allowed-rules',
                'final_filter_snps', 'final_filter_indels'], cwd=ROOT, env=env,
                capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            for kind, expected in (('snps', 5000), ('indels', 0)):
                path = out/'final'/f'filtered_{kind}.final.vcf.gz'
                with gzip.open(path, 'rt') as f:
                    self.assertEqual(sum(1 for line in f if not line.startswith('#')), expected)
                self.assertGreater(Path(str(path)+'.tbi').stat().st_size, 0)
