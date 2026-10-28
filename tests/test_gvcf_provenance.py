"""Regression tests for reference validation and safe GVCF adoption."""
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

SCRIPT = Path(__file__).resolve().parents[1] / 'workflow/scripts/validate_gvcf_provenance.py'
spec = importlib.util.spec_from_file_location('provenance', SCRIPT)
provenance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ref = self.root / 'reference.fa'
        self.ref.write_text('>chr1\nACGT\n')
        self.fai = self.root / 'reference.fa.fai'
        self.fai.write_text('chr1\t4\t6\t4\t5\n')
        self.sha = self.root / 'reference.sha256'
        self.sha.write_text(provenance.sha256_file(self.ref))
        self.gvcf = self.root / 'sample.g.vcf.gz'
        self.meta = self.root / 'sample.g.vcf.gz.meta.json'
        self.job = SimpleNamespace(
            input=SimpleNamespace(gvcf=self.gvcf, ref=self.ref, fai=self.fai, reference_sha256=self.sha),
            output=SimpleNamespace(metadata=self.meta), log=[self.root / 'validation.log'],
            params=SimpleNamespace(allow_header_only='false'),
            wildcards=SimpleNamespace(sample='sample'))

    def write_gvcf(self, sample='sample', reference=None, contig_length=4, gatk=False):
        reference = reference or self.ref
        header = (f'##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --reference {reference}",Version="4">\n'
                  if gatk else f'##reference=file://{reference}\n')
        with gzip.open(self.gvcf, 'wt') as f:
            f.write('##fileformat=VCFv4.2\n' + header)
            f.write(f'##contig=<ID=chr1,length={contig_length}>\n')
            f.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n')
            f.write('chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=4\tGT\t0/0\n')

    def test_gatk_command_reference_is_accepted(self):
        self.write_gvcf(gatk=True)
        provenance.validate(self.job)
        meta = json.loads(self.meta.read_text())
        self.assertEqual(meta['reference_sha256'], self.sha.read_text())
        self.assertEqual(meta['gvcf_sha256'], provenance.sha256_file(self.gvcf))

    def test_uncompressed_legacy_gvcf_is_accepted(self):
        self.write_gvcf()
        plain = self.root / 'sample.snps.indels.g.vcf'
        with gzip.open(self.gvcf, 'rb') as source:
            plain.write_bytes(source.read())
        self.job.input.gvcf = plain
        provenance.validate(self.job)
        self.assertEqual(json.loads(self.meta.read_text())['gvcf'], str(plain))

    def test_wrong_sample_is_rejected(self):
        self.write_gvcf(sample='other')
        with self.assertRaisesRegex(ValueError, 'sample mismatch'):
            provenance.validate(self.job)

    def test_same_lengths_different_sequence_is_rejected(self):
        other = self.root / 'other.fa'
        other.write_text('>chr1\nAAAA\n')
        self.write_gvcf(reference=other)
        with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
            provenance.validate(self.job)

    def test_wrong_contigs_are_rejected(self):
        self.write_gvcf(contig_length=5)
        with self.assertRaisesRegex(ValueError, 'contig mismatch'):
            provenance.validate(self.job)

    def test_truncated_gzip_is_rejected(self):
        self.write_gvcf()
        self.gvcf.write_bytes(self.gvcf.read_bytes()[:-8])
        with self.assertRaises(EOFError):
            provenance.validate(self.job)

    def test_sidecar_remains_valid_after_legacy_fallback_disabled(self):
        self.write_gvcf(reference=self.root / 'missing.fa')
        self.job.params.allow_header_only = 'true'
        provenance.validate(self.job)
        self.job.params.allow_header_only = 'false'
        self.job.params.previous_metadata = json.loads(self.meta.read_text())
        self.meta.unlink()  # model Snakemake removing an output before a rerun
        provenance.validate(self.job)
        self.assertEqual(json.loads(self.meta.read_text())['validation_mode'], 'validated_sidecar')
        self.write_gvcf(sample='sample', reference=self.root / 'missing.fa', gatk=True)
        with self.assertRaisesRegex(ValueError, 'Cannot prove'):
            provenance.validate(self.job)


if __name__ == '__main__':
    unittest.main()
