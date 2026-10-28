"""Run both population-analysis branches on a synthetic VCF using installed tools."""
import csv
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = os.environ.get('SNAKEMAKE_BINARY') or shutil.which('snakemake')
ANALYSIS_BIN = os.environ.get('POSTPROCESSING_TEST_BIN')


@unittest.skipUnless(SNAKEMAKE and ANALYSIS_BIN, 'Set POSTPROCESSING_TEST_BIN for PLINK/R/LEA integration')
class PopulationIntegrationTests(unittest.TestCase):
    def test_figures_lists_and_model_reuse(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            final = root/'results/final'
            final.mkdir(parents=True)
            samples = [f'sample{i:02d}' for i in range(60)]
            plain = root/'input.vcf'
            rng = random.Random(7634)
            with plain.open('w') as f:
                f.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n'
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n')
                for i in range(600):
                    genotypes = [rng.choice(['0/0','0/1','1/1']) for _ in samples]
                    f.write(f'1\t{i+1}\t.\tA\tG\t60\tPASS\t.\tGT\t'+'\t'.join(genotypes)+'\n')
            env = dict(os.environ, PATH=ANALYSIS_BIN+os.pathsep+os.environ['PATH'],
                       XDG_CACHE_HOME=str(root/'cache'))
            compressed = final/'filtered_snps.final.vcf.gz'
            subprocess.run([str(Path(ANALYSIS_BIN)/'bcftools'),'view','-Oz','-o',str(compressed),str(plain)],
                           check=True, env=env, capture_output=True)
            subprocess.run([str(Path(ANALYSIS_BIN)/'bcftools'),'index','-t',str(compressed)], check=True, env=env)
            clades, pca_list, tree_list = root/'clades.tsv', root/'pca.txt', root/'tree.txt'
            clades.write_text('sample\tclade\n'+''.join(f'{sample}\tgroup{i%3}\n' for i,sample in enumerate(samples)))
            pca_list.write_text('sample00\nsample01\n')
            tree_list.write_text('sample02\nsample03\nsample04\n')
            settings = {
                'enabled': True, 'prefix': 'synthetic', 'metadata': {'clades': str(clades)},
                'pca': {'components':3, 'sample_list':str(pca_list)},
                'phylogenetic_tree': {'layout':'rectangular','width':4,'height':4,'dpi':72,'sample_list':str(tree_list)},
                'lea': {'k_min':2,'k_max':3,'plot_k_min':2,'plot_k_max':3,'repetitions':1,
                        'width':4,'height':3,'multi_height':4},
            }
            config = root/'config.json'
            conf = {'ref':'tests/minimal/reference.fa','samples_tsv':'tests/minimal/samples.tsv',
                    'outdir':str(root/'results'), 'resources':{'postprocessing':{'threads':1}}, 'postprocessing':settings}
            command = [SNAKEMAKE,'postprocessing','--snakefile','workflow/Snakefile','--configfile',str(config),
                       '--cores','2','--allowed-rules','postprocessing','postprocess_pca_tree_prepare',
                       'postprocess_pca_tree_plot','postprocess_lea_prepare','postprocess_lea_model','postprocess_lea_plot']
            def execute():
                config.write_text(json.dumps(conf))
                result = subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,timeout=180)
                if result.returncode:
                    logs = '\n'.join(p.read_text(errors='replace')[-3000:]
                                     for p in (root/'results/logs/postprocessing').glob('*.log'))
                    self.fail(result.stdout+result.stderr+'\n'+logs)
                return result.stdout+result.stderr
            execute()
            tree = root/'results/postprocessing/pca_tree'
            model = root/'results/postprocessing/lea/synthetic.LEA.snmfProject'
            before = model.stat().st_mtime_ns
            with (tree/'synthetic.tree.pca_scores.tsv').open() as f:
                self.assertEqual({r['IID'] for r in csv.DictReader(f,delimiter='\t')},{'sample00','sample01'})
            info = dict(line.split('\t',1) for line in (tree/'synthetic.tree.tree_info.tsv').read_text().splitlines())
            self.assertEqual(info['tips'],'3')
            settings['lea']['best_k'] = 3
            execute()
            self.assertEqual(model.stat().st_mtime_ns,before)
            eigenvec = tree/'synthetic.pca.eigenvec'
            eigen_time = eigenvec.stat().st_mtime_ns
            pca_list.write_text('sample05\n')
            execute()
            self.assertEqual(eigenvec.stat().st_mtime_ns,eigen_time)
            self.assertEqual(model.stat().st_mtime_ns,before)
            settings['lea']['seed'] = 9876
            execute()
            self.assertNotEqual(model.stat().st_mtime_ns,before)
            self.assertIn('Nothing to be done',execute())


if __name__ == '__main__':
    unittest.main()
