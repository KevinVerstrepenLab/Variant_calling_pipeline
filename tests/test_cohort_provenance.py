"""Regression coverage for parallel cohort validation."""
import gzip
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('cohort', ROOT / 'workflow/scripts/validate_cohort_provenance.py')
COHORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COHORT)


class CohortTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        ref = self.root / 'reference.fa'
        ref.write_text('>chr1\nACGT\n')
        fai = self.root / 'reference.fa.fai'
        fai.write_text('chr1\t4\t6\t4\t5\n')
        sha = self.root / 'reference.sha256'
        sha.write_text(hashlib.sha256(ref.read_bytes()).hexdigest())
        paths = []
        for sample, suffix in [('one', '.snps.indels.g.vcf'), ('two', '.g.vcf.gz')]:
            path = self.root / (sample + suffix)
            opener = gzip.open if suffix.endswith('.gz') else open
            with opener(path, 'wt') as f:
                f.write(f'##fileformat=VCFv4.2\n##reference={ref}\n')
                f.write('##contig=<ID=chr1,length=4>\n')
                f.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n')
            paths.append(str(path))
        self.job = SimpleNamespace(
            input=SimpleNamespace(validator=str(ROOT / 'workflow/scripts/validate_gvcf_provenance.py'),
                                  gvcfs=paths, ref=str(ref), fai=str(fai), reference_sha256=str(sha)),
            output=SimpleNamespace(manifest=str(self.root / 'cohort.validated.json')),
            params=SimpleNamespace(samples=['one', 'two'], previous_metadata={},
                                   metadata=[str(self.root / f'{s}.meta.json') for s in ['one', 'two']],
                                   logs_dir=str(self.root / 'logs'), allow_header_only='false'),
            threads=2, log=[str(self.root / 'cohort.log')])

    def test_all_formats_receive_their_own_metadata(self):
        COHORT.validate_cohort(self.job)
        for sample, output in zip(self.job.params.samples, self.job.params.metadata):
            self.assertEqual(json.loads(Path(output).read_text())['sample'], sample)
        manifest = json.loads(Path(self.job.output.manifest).read_text())
        self.assertEqual(len(manifest['samples']), 2)
        self.assertIn('Validated 2/2 GVCFs', Path(self.job.log[0]).read_text())

    def test_one_failed_sample_fails_the_cohort(self):
        self.job.params.samples[0] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'failed for wrong'):
            COHORT.validate_cohort(self.job)
        self.assertIn('FAILED wrong', Path(self.job.log[0]).read_text())


if __name__ == '__main__':
    unittest.main()
