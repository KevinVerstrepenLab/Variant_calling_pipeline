"""Coverage screening: synthetic gains/losses and an optional real BAM-to-figure run."""
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow'))
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from configuration import coverage_settings
from coverage_windows import write_windows
from plot_coverage import (analyze, gc_correct, smooth_contiguous, weighted_median,
                           display_layout, display_limits, write_results)
import numpy as np

SNAKEMAKE = os.environ.get('SNAKEMAKE_BINARY') or shutil.which('snakemake')
COVERAGE_BIN = os.environ.get('COVERAGE_TEST_BIN')


class CoverageTests(unittest.TestCase):
    def test_display_defaults_and_reference_order_wrapping(self):
        settings = coverage_settings({'coverage_plots': {'plot_min_contig_bp': None}})
        rows = [{'chromosome': f'chr{i}', 'start': '0', 'end': '50000'} for i in range(17)]
        rows.insert(3, {'chromosome': 'short', 'start': '0', 'end': '49999'})
        bands, lengths = display_layout(rows, settings)
        self.assertEqual(settings['plot_min_contig_bp'], 50000)
        self.assertEqual([len(band) for band in bands], [16, 1])
        self.assertEqual(bands[-1], ['chr16'])
        self.assertNotIn('short', [name for band in bands for name in band])
        self.assertEqual(sum(map(len, display_layout(rows, dict(settings, plot_min_contig_bp=0))[0])), 18)
        with self.assertRaisesRegex(ValueError, 'No sequences meet'):
            display_layout(rows, dict(settings, plot_min_contig_bp=100000))

    def test_isolated_extremes_do_not_set_plot_scale(self):
        depth = np.array([0.] + [20.]*100 + [1e9])
        limits = display_limits(depth, [0.01, 0.99], anchors=[20], minimum_span=4, nonnegative=True)
        self.assertGreater(limits[0], 0)
        self.assertLess(limits[1], 30)
        self.assertGreater(display_limits(depth, None)[1], 1e9)

    def test_wrapped_plot_keeps_hidden_sequences_and_outliers_in_tables(self):
        with tempfile.TemporaryDirectory() as temp:
            settings = coverage_settings({'coverage_plots': {'width': 8, 'height': 4, 'dpi': 50}})
            rows = []
            for i in range(18):
                length = 40000 if i == 17 else 50000
                for start in range(0, length, 10000):
                    rows.append({'chromosome': f'chr{i}', 'start': str(start), 'end': str(start+10000),
                                 'gc': '0.5', 'n_fraction': '0', 'eligible': '1'})
            depth = np.full(len(rows), 20.)
            depth[0], depth[1] = 1e8, 0
            ratios = np.log2((depth+0.5)/20.5)
            result = (rows, depth, depth.copy(), ratios, ratios.copy(), np.zeros(len(rows)), [],
                      {'baseline_depth': 20, 'low_depth': False, 'gc_correction': 'disabled'})
            prefix = Path(temp)/'plot'
            write_results(prefix, 'synthetic', result, settings)
            qc = json.loads(Path(str(prefix)+'.qc.json').read_text())['plot']
            self.assertEqual(qc['sequence_rows'], 2)
            self.assertEqual(qc['hidden_contigs'], ['chr17'])
            self.assertLess(qc['depth_ylim'][1], 30)
            self.assertEqual(qc['rows'][0]['depth_outside_range'], {'below': 1, 'above': 1})
            with gzip.open(str(prefix)+'.windows.tsv.gz', 'rt') as f:
                saved = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(len(saved), len(rows))
            self.assertEqual(float(saved[0]['mean_depth']), 1e8)
            self.assertEqual(saved[-1]['chromosome'], 'chr17')
            for suffix in ('.png', '.pdf'):
                self.assertGreater(Path(str(prefix)+suffix).stat().st_size, 0)

    def test_coordinates_and_ambiguous_mask(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'ref.fa').write_text('>chr1\nACGTNNACGT\n>tiny\nA\n')
            settings = coverage_settings({'coverage_plots': {'window_bp': 4, 'min_contig_bp': 2}})
            write_windows(root/'ref.fa', root/'windows.bed', root/'windows.tsv', settings)
            self.assertEqual((root/'windows.bed').read_text(), 'chr1\t0\t4\nchr1\t4\t8\nchr1\t8\t10\n')
            with (root/'windows.tsv').open() as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual([r['eligible'] for r in rows], ['1', '0', '1'])
            self.assertEqual(float(rows[0]['gc']), 0.5)

    def test_local_gains_losses_and_zero_depth(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'ref.fa').write_text('>chr1\n'+'A'*30000+'\n')
            settings = coverage_settings({'coverage_plots': {'window_bp': 1000, 'gc_correct': False, 'smooth_windows': 3}})
            write_windows(root/'ref.fa', root/'windows.bed', root/'windows.tsv', settings)
            depths = [10]*10 + [20]*5 + [0]*5 + [10]*10
            with gzip.open(root/'regions.gz', 'wt') as f:
                for i, depth in enumerate(depths):
                    f.write(f'chr1\t{i*1000}\t{(i+1)*1000}\t{depth}\n')
            result = analyze(root/'regions.gz', root/'windows.tsv', settings)
            self.assertEqual(result[-1]['baseline_depth'], 10)
            self.assertEqual({r[3] for r in result[-2]}, {'gain', 'loss'})
            self.assertTrue(np.isfinite(result[4]).all())
            self.assertLess(result[4][16], -4)

    def test_weighting_and_contiguous_smoothing(self):
        self.assertEqual(weighted_median([2,10,20], [1,10,1]), 10)
        values = np.array([1.,np.nan,-1.,3.])
        result = smooth_contiguous(values, ['a','a','a','b'], [0,10,20,0], [10,20,30,10], [True,False,True,True], 5)
        np.testing.assert_allclose(result, values, equal_nan=True)

    def test_optional_gc_correction(self):
        gc = np.linspace(0.2,0.8,1000)
        depth = 10+30*gc
        settings = coverage_settings({})
        corrected, status = gc_correct(depth,gc,np.ones(1000),np.ones(1000,dtype=bool),settings)
        self.assertIn('applied', status)
        self.assertLess(np.std(corrected), np.std(depth)/10)
        raw, status = gc_correct(depth,gc,np.ones(1000),np.ones(1000,dtype=bool),dict(settings,gc_correct=False))
        np.testing.assert_array_equal(raw, depth)

    def test_invalid_parameters(self):
        for override in ({'window_bp':0}, {'smooth_windows':2}, {'pseudocount':0}, {'max_n_fraction':2},
                         {'plot_min_contig_bp': -1}, {'chromosomes_per_row': 0},
                         {'outlier_quantiles': [0.99, 0.01]}, {'outlier_quantiles': [0, 2]}):
            with self.assertRaises(ValueError): coverage_settings({'coverage_plots':override})


@unittest.skipUnless(SNAKEMAKE, 'Snakemake required')
class CoverageDagTests(unittest.TestCase):
    def test_existing_bams_and_new_samples_without_recalling_reused_gvcfs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'bam').mkdir(); (root/'gvcf').mkdir()
            (root/'bam/old.sorted.bam').write_bytes(b'placeholder for dry run')
            with gzip.open(root/'gvcf/paired.g.vcf.gz','wt') as f:
                f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpaired\n')
            config=root/'override.json'
            config.write_text(json.dumps({'outdir':str(root),'coverage_plots':{'enabled':True}}))
            result=subprocess.run([SNAKEMAKE,'coverage_plots','--snakefile','workflow/Snakefile','--configfile',
                                   'tests/minimal/config.yaml',str(config),'--cores','2','--dry-run','--printshellcmds'],cwd=ROOT,
                                  capture_output=True,text=True,env=dict(os.environ,XDG_CACHE_HOME=str(root/'cache')))
            output=result.stdout+result.stderr
            self.assertEqual(result.returncode,0,output)
            self.assertIn('old.sorted.bam',output)
            self.assertIn('single.sorted.bam',output)
            self.assertNotIn('paired.sorted.bam',output)
            self.assertNotIn('haplotype_caller_gvcf',output)
            self.assertNotIn('final_filter_snps',output)


@unittest.skipUnless(SNAKEMAKE and COVERAGE_BIN, 'Set COVERAGE_TEST_BIN for real mosdepth/samtools testing')
class CoverageIntegrationTests(unittest.TestCase):
    def test_bam_to_figures_and_reuse(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); bams=root/'input bams'; bams.mkdir()
            reference=root/'ref.fa'; reference.write_text('>chr1\n'+'ACGT'*2500+'\n')
            sam=root/'reads.sam'
            with sam.open('w') as f:
                f.write('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n')
                for start in range(0,10000,100):
                    copies=10 if 2000<=start<4000 else 30 if 6000<=start<8000 else 0 if start>=9000 else 20
                    for copy in range(copies):
                        f.write(f'r{start}_{copy}\t0\tchr1\t{start+1}\t60\t100M\t*\t0\t0\t'+'ACGT'*25+'\t'+'I'*100+'\n')
            env=dict(os.environ,PATH=str(Path(sys.executable).parent)+os.pathsep+COVERAGE_BIN+os.pathsep+os.environ['PATH'],
                     XDG_CACHE_HOME=str(root/'cache'),MPLCONFIGDIR=str(root/'mpl'))
            bam=bams/'synthetic.sorted.bam'
            subprocess.run([str(Path(COVERAGE_BIN)/'samtools'),'view','-b','-o',str(bam),str(sam)],check=True,env=env)
            out=root/'results'
            settings={'enabled':True,'bam_dir':str(bams),'window_bp':1000,'min_contig_bp':1,'gc_correct':False,
                      'smooth_windows':1,'min_region_windows':1,'width':8,'height':5,'dpi':72,
                      'plot_min_contig_bp':0}
            config=root/'config.json'
            conf={'ref':str(reference),'samples_tsv':'tests/minimal/samples.tsv','outdir':str(out),'coverage_plots':settings}
            config.write_text(json.dumps(conf))
            command=[SNAKEMAKE,'coverage_plots','--snakefile','workflow/Snakefile','--configfile',str(config),'--cores','2']
            def execute():
                result=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,timeout=120)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                return result.stdout+result.stderr
            execute()
            prefix=out/'CoveragePlots/synthetic'
            for suffix in ('.png','.pdf','.windows.tsv.gz','.candidate_regions.tsv','.qc.json'):
                self.assertGreater(Path(str(prefix)+suffix).stat().st_size,0)
            self.assertEqual(json.loads(Path(str(prefix)+'.qc.json').read_text())['baseline_depth'],20)
            with gzip.open(str(prefix)+'.windows.tsv.gz','rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),10)
            self.assertEqual(float(rows[-1]['mean_depth']),0)
            self.assertLess(float(rows[-1]['log2_ratio']),-5)
            with open(str(prefix)+'.candidate_regions.tsv') as f:
                self.assertEqual({r['direction'] for r in csv.DictReader(f,delimiter='\t')},{'gain','loss'})
            depth=out/'CoveragePlots/depth/synthetic.regions.bed.gz'; before=depth.stat().st_mtime_ns
            settings['log2_threshold']=0.5; config.write_text(json.dumps(conf)); execute()
            self.assertEqual(depth.stat().st_mtime_ns,before)
            settings['outlier_quantiles'] = None
            settings['chromosomes_per_row'] = 8
            settings['plot_min_contig_bp'] = 5000
            config.write_text(json.dumps(conf)); execute()
            self.assertEqual(depth.stat().st_mtime_ns, before)
            self.assertIn('Nothing to be done',execute())
            old_time = bam.stat().st_mtime_ns
            os.utime(bam, ns=(old_time+1000000000, old_time+1000000000))
            self.assertIn('coverage_depth', execute())
            self.assertNotEqual(depth.stat().st_mtime_ns, before)
            self.assertFalse(Path(str(bam)+'.bai').exists())


if __name__=='__main__': unittest.main()
