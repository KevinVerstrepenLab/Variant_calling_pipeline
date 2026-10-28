import json

POSTPROCESS_ENGINE = WORKFLOW / "scripts" / "postprocess_vcf.py"
POSTPROCESS_RUNNER = WORKFLOW / "scripts" / "run_postprocessing.py"
POSTPROCESSING_TARGETS = [path for stage in ("pca_tree", "lea")
                          for path in postprocessing_outputs(POSTPROCESSING, OUT, stage)]


rule postprocessing:
    input:
        POSTPROCESSING_TARGETS


if POSTPROCESSING["enabled"] and POSTPROCESSING["pca_tree"]["enabled"]:

    rule postprocess_pca_tree_prepare:
        input:
            vcf=f"{FINAL}/filtered_snps.final.vcf.gz",
            index=f"{FINAL}/filtered_snps.final.vcf.gz.tbi",
            engine=str(POSTPROCESS_ENGINE),
            runner=str(POSTPROCESS_RUNNER),
            metadata=[POSTPROCESSING["pca_tree"]["keep"]] if POSTPROCESSING["pca_tree"].get("keep") else [],
        output:
            postprocessing_outputs(POSTPROCESSING, OUT, "pca_tree", "prepare")
        params:
            settings=json.dumps(postprocessing_phase_settings(POSTPROCESSING, "pca_tree", "prepare"), sort_keys=True),
            prefix=POSTPROCESSING["prefix"],
            outdir=lambda w, output: str(Path(output[0]).parent)
        threads: config.get("postprocessing_threads", 8)
        resources:
            mem_mb=config.get("postprocessing_mem_mb", 16000)
        conda: str(ENVS / "postprocessing.yaml")
        log: f"{LOGS}/postprocessing/pca_tree_prepare.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {log:q})"
            python {input.runner:q} --stage pca_tree --phase prepare --input {input.vcf:q} \
                --settings {params.settings:q} --prefix {params.prefix:q} \
                --outdir {params.outdir:q} --threads {threads} --mem-mb {resources.mem_mb} > {log:q} 2>&1
            """

    rule postprocess_pca_tree_plot:
        input:
            vcf=f"{FINAL}/filtered_snps.final.vcf.gz",
            index=f"{FINAL}/filtered_snps.final.vcf.gz.tbi",
            engine=str(POSTPROCESS_ENGINE),
            runner=str(POSTPROCESS_RUNNER),
            prepare=postprocessing_outputs(POSTPROCESSING, OUT, "pca_tree", "prepare"),
            metadata=[POSTPROCESSING["pca_tree"][key] for key in ("keep", "rename", "clades", "clade_colors", "pca_sample_list", "tree_sample_list") if POSTPROCESSING["pca_tree"].get(key)],
            rscript=str(WORKFLOW / "scripts" / "build_tree_postprocess.R"),
            palette=str(WORKFLOW / "scripts" / "clade_palette_utils.R"),
        output:
            postprocessing_outputs(POSTPROCESSING, OUT, "pca_tree", "plot")
        params:
            settings=json.dumps(postprocessing_phase_settings(POSTPROCESSING, "pca_tree", "plot"), sort_keys=True),
            prefix=POSTPROCESSING["prefix"],
            outdir=lambda w, output: str(Path(output[0]).parent)
        threads: config.get("postprocessing_threads", 8)
        resources:
            mem_mb=config.get("postprocessing_mem_mb", 16000)
        conda: str(ENVS / "postprocessing.yaml")
        log: f"{LOGS}/postprocessing/pca_tree_plot.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {log:q})"
            python {input.runner:q} --stage pca_tree --phase plot --input {input.vcf:q} \
                --settings {params.settings:q} --prefix {params.prefix:q} \
                --outdir {params.outdir:q} --threads {threads} --mem-mb {resources.mem_mb} > {log:q} 2>&1
            """


if POSTPROCESSING["enabled"] and POSTPROCESSING["lea"]["enabled"]:

    rule postprocess_lea_prepare:
        input:
            vcf=f"{FINAL}/filtered_snps.final.vcf.gz",
            index=f"{FINAL}/filtered_snps.final.vcf.gz.tbi",
            engine=str(POSTPROCESS_ENGINE),
            runner=str(POSTPROCESS_RUNNER),
        output:
            postprocessing_outputs(POSTPROCESSING, OUT, "lea", "prepare")
        params:
            settings=json.dumps(postprocessing_phase_settings(POSTPROCESSING, "lea", "prepare"), sort_keys=True),
            prefix=POSTPROCESSING["prefix"],
            outdir=lambda w, output: str(Path(output[0]).parent)
        threads: config.get("postprocessing_threads", 8)
        resources:
            mem_mb=config.get("postprocessing_mem_mb", 16000)
        conda: str(ENVS / "postprocessing.yaml")
        log: f"{LOGS}/postprocessing/lea_prepare.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {log:q})"
            python {input.runner:q} --stage lea --phase prepare --input {input.vcf:q} \
                --settings {params.settings:q} --prefix {params.prefix:q} \
                --outdir {params.outdir:q} --threads {threads} --mem-mb {resources.mem_mb} > {log:q} 2>&1
            """
    if POSTPROCESSING["lea"].get("run_snmf", True):

        rule postprocess_lea_model:
            input:
                vcf=f"{FINAL}/filtered_snps.final.vcf.gz",
                index=f"{FINAL}/filtered_snps.final.vcf.gz.tbi",
                engine=str(POSTPROCESS_ENGINE),
                runner=str(POSTPROCESS_RUNNER),
                prepare=postprocessing_outputs(POSTPROCESSING, OUT, "lea", "prepare"),
                rscript=str(WORKFLOW / "scripts" / "run_lea.R"),
            output:
                postprocessing_outputs(POSTPROCESSING, OUT, "lea", "model")
            params:
                settings=json.dumps(postprocessing_phase_settings(POSTPROCESSING, "lea", "model"), sort_keys=True),
                prefix=POSTPROCESSING["prefix"],
                outdir=lambda w, output: str(Path(output[0]).parent)
            threads: config.get("postprocessing_threads", 8)
            resources:
                mem_mb=config.get("postprocessing_mem_mb", 16000)
            conda: str(ENVS / "postprocessing.yaml")
            log: f"{LOGS}/postprocessing/lea_model.log"
            shell:
                """
                set -euo pipefail
                mkdir -p "$(dirname {log:q})"
                python {input.runner:q} --stage lea --phase model --input {input.vcf:q} \
                    --settings {params.settings:q} --prefix {params.prefix:q} \
                    --outdir {params.outdir:q} --threads {threads} --mem-mb {resources.mem_mb} > {log:q} 2>&1
                """
    if POSTPROCESSING["lea"].get("run_snmf", True):

        rule postprocess_lea_plot:
            input:
                vcf=f"{FINAL}/filtered_snps.final.vcf.gz",
                index=f"{FINAL}/filtered_snps.final.vcf.gz.tbi",
                engine=str(POSTPROCESS_ENGINE),
                runner=str(POSTPROCESS_RUNNER),
                prepare=postprocessing_outputs(POSTPROCESSING, OUT, "lea", "prepare"),
                model=postprocessing_outputs(POSTPROCESSING, OUT, "lea", "model"),
                metadata=[POSTPROCESSING["lea"][key] for key in ("keep", "rename", "clades", "clade_colors") if POSTPROCESSING["lea"].get(key)],
                rscript=str(WORKFLOW / "scripts" / "run_lea.R"),
                palette=str(WORKFLOW / "scripts" / "clade_palette_utils.R"),
            output:
                postprocessing_outputs(POSTPROCESSING, OUT, "lea", "plot")
            params:
                settings=json.dumps(postprocessing_phase_settings(POSTPROCESSING, "lea", "plot"), sort_keys=True),
                prefix=POSTPROCESSING["prefix"],
                outdir=lambda w, output: str(Path(output[0]).parent)
            threads: config.get("postprocessing_threads", 8)
            resources:
                mem_mb=config.get("postprocessing_mem_mb", 16000)
            conda: str(ENVS / "postprocessing.yaml")
            log: f"{LOGS}/postprocessing/lea_plot.log"
            shell:
                """
                set -euo pipefail
                mkdir -p "$(dirname {log:q})"
                python {input.runner:q} --stage lea --phase plot --input {input.vcf:q} \
                    --settings {params.settings:q} --prefix {params.prefix:q} \
                    --outdir {params.outdir:q} --threads {threads} --mem-mb {resources.mem_mb} > {log:q} 2>&1
                """
