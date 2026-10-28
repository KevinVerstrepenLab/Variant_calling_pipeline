"""Render synthetic PCA/tree subsets with an explicitly selected analysis R environment."""
import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get('POSTPROCESSING_RSCRIPT')


@unittest.skipUnless(RSCRIPT, 'Set POSTPROCESSING_RSCRIPT to an Rscript with plotting dependencies')
class SampleListPlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = self.root / 'figure'
        samples = [f'sample{i}' for i in range(1, 7)]
        (self.root / 'ids').write_text(''.join(f'{sample}\t{sample}\n' for sample in samples))
        matrix = [[0 if i == j else abs(i-j)+1 for j in range(6)] for i in range(6)]
        (self.root / 'dist').write_text(''.join('\t'.join(map(str, row))+'\n' for row in matrix))
        (self.root / 'eigenvec').write_text('#FID\tIID\tPC1\tPC2\n'+''.join(
            f'{sample}\t{sample}\t{i}\t{i*i}\n' for i, sample in enumerate(samples, 1)))
        (self.root / 'eigenval').write_text('10\n5\n')
        (self.root / 'rename').write_text('sample1\tRenamed1\n')
        self.command = [RSCRIPT, str(ROOT / 'workflow/scripts/build_tree_postprocess.R'),
                        '--dist', str(self.root / 'dist'), '--ids', str(self.root / 'ids'),
                        '--pca_eigenvec', str(self.root / 'eigenvec'),
                        '--pca_eigenval', str(self.root / 'eigenval'),
                        '--rename', str(self.root / 'rename'), '--out', str(self.out),
                        '--layout', 'rectangular', '--width', '4', '--height', '4', '--dpi', '72']

    def render(self, pca=None, tree=None):
        command = list(self.command)
        for method, samples in [('pca', pca), ('tree', tree)]:
            if samples is not None:
                path = self.root / f'{method}.txt'
                path.write_text(samples)
                command += [f'--{method}_sample_list', str(path)]
        return subprocess.run(command, capture_output=True, text=True, timeout=90)

    def test_independent_subsets_and_renamed_ids(self):
        result = self.render(pca='# PCA selection\nRenamed1\nsample3\nsample3\n\n',
                             tree='sample2\nsample4\nsample5\n')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with Path(str(self.out)+'.pca_scores.tsv').open() as handle:
            scores = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([row['IID'] for row in scores], ['sample1', 'sample3'])
        self.assertEqual([float(row['PC1']) for row in scores], [1, 3])
        info = dict(line.rstrip('\n').split('\t', 1) for line in Path(str(self.out)+'.tree_info.tsv').read_text().splitlines())
        self.assertEqual(info['tips'], '3')
        self.assertEqual(info['inferred_tips'], '6')
        self.assertEqual(set(info['displayed_samples'].split(',')), {'sample2', 'sample4', 'sample5'})
        self.assertIn('Renamed1', Path(str(self.out)+'.nwk').read_text())
        for suffix in ('.pca.pdf', '.pca.png', '-rectangular.pdf', '-rectangular.png'):
            self.assertGreater(Path(str(self.out)+suffix).stat().st_size, 0)

    def test_nonfinite_entropy_preserves_snmf_run_numbers(self):
        expression = r'''expr <- parse(file=commandArgs(TRUE)[1]);
for (e in expr) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.symbol(e[[2]]) &&
      as.character(e[[2]]) == "cross_entropy_values") eval(e)
}
v <- cross_entropy_values(c(NA_real_, 0.2, 0.1, Inf));
stopifnot(length(v) == 4, which.min(v) == 3, is.na(v[1]), is.na(v[4]))'''
        result = subprocess.run([RSCRIPT, '-e', expression, str(ROOT / 'workflow/scripts/run_lea.R')],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_no_matching_pca_samples_fails(self):
        result = self.render(pca='not_in_the_cohort\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('PCA sample_list matched no samples', result.stderr)

    def test_tree_requires_two_matching_tips(self):
        result = self.render(tree='sample1\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('TREE sample_list requires at least two matching samples', result.stderr)


if __name__ == '__main__':
    unittest.main()
