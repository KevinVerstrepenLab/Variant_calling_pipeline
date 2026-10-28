"""Exercise public GATK filters with scientific notation and missing INFO annotations."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
GATK = os.environ.get('GATK_TEST_BINARY') or shutil.which('gatk')


@unittest.skipUnless(GATK, 'GATK required')
class HardFilterTests(unittest.TestCase):
    def test_scientific_notation_and_missing_annotations(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ref = root/'ref.fa'
            ref.write_text('>chr1\n'+'A'*200+'\n')
            (root/'ref.fa.fai').write_text('chr1\t200\t6\t200\t201\n')
            (root/'ref.dict').write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:200\n')
            cfg = yaml.safe_load((ROOT/'config/config.yaml').read_text())
            keys = ('QD', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'SOR')
            header = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n'
            header += ''.join(f'##INFO=<ID={key},Number=1,Type=Float,Description="Test">\n' for key in keys)
            header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            inputs = {
                'snp': ['QD=10;FS=1;MQ=60;MQRankSum=-1.593e+00;ReadPosRankSum=0;SOR=1',
                        'QD=2', 'MQRankSum=-9.0', '.'],
                'indel': ['QD=10;FS=1;ReadPosRankSum=-1.593e+00;SOR=1',
                          'QD=2', 'ReadPosRankSum=-11.0', '.'],
            }
            env = dict(os.environ, PATH=str(Path(GATK).absolute().parent)+os.pathsep+os.environ['PATH'])
            for kind, annotations in inputs.items():
                with self.subTest(kind=kind):
                    vcf, out = root/f'{kind}.vcf', root/f'{kind}.filtered.vcf'
                    alt = 'C' if kind == 'snp' else 'AC'
                    vcf.write_text(header+''.join(f'chr1\t{i}\t.\tA\t{alt}\t60\tPASS\t{info}\n'
                                                for i, info in enumerate(annotations, 1)))
                    result = subprocess.run([GATK, '--java-options', '-Xms128m -Xmx512m', 'VariantFiltration',
                        '-R', str(ref), '-V', str(vcf), '--filter-expression', cfg['filtering']['hard'][kind],
                        '--filter-name', 'test_filter', '-O', str(out)], env=env, capture_output=True, text=True, timeout=90)
                    self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                    filters = [line.split('\t')[6] for line in out.read_text().splitlines() if not line.startswith('#')]
                    self.assertEqual(filters, ['PASS', 'test_filter', 'test_filter', 'PASS'])
