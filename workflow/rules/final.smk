# Final SNP/INDEL filtering with AD-masking, sample/site missingness.
# Requires variables from Snakefile: ENVS, LOGS, FINAL, SNPS_DIR, INDELS_DIR.

def _ad_expr():
    ad = config.get("ad_filter", {})
    min_dp = int(ad.get("min_dp", 5))
    lo = float(ad.get("ambig_low", 0.1))
    hi = float(ad.get("ambig_high", 0.9))
    # For biallelic records: mask if DP < min_dp OR ALT fraction in (lo, hi)
    return (
        f'FMT/DP<{min_dp} || '
        f'((FMT/AD[0]+FMT/AD[1])>0 && '
        f'(FMT/AD[1])/(FMT/AD[0]+FMT/AD[1])>{lo} && '
        f'(FMT/AD[1])/(FMT/AD[0]+FMT/AD[1])<{hi})'
    )

# ---------------- SNPs ----------------
rule final_filter_snps:
    input:
        vcf = f"{SNPS_DIR}/filtered_snps.merged.vcf.gz",
        tbi = f"{SNPS_DIR}/filtered_snps.merged.vcf.gz.tbi"
    output:
        passvcf = temp(f"{FINAL}/snps.step1_filtered.vcf.gz"),
        passtbi = temp(f"{FINAL}/snps.step1_filtered.vcf.gz.tbi"),
        advcf   = temp(f"{FINAL}/snps.step2_adclean.vcf.gz"),
        adtbi   = temp(f"{FINAL}/snps.step2_adclean.vcf.gz.tbi"),
        step3vcf = temp(f"{FINAL}/snps.step3_filtered.vcf.gz"),
        step3tbi = temp(f"{FINAL}/snps.step3_filtered.vcf.gz.tbi"),
        removes = temp(f"{FINAL}/snps.remove_samples.txt"),
        vcf = f"{FINAL}/filtered_snps.final.vcf.gz",
        tbi = f"{FINAL}/filtered_snps.final.vcf.gz.tbi"
    conda: f"{ENVS}/vcfops.yaml"
    threads: 4
    log: f"{LOGS}/final/snps.log"
    params:
        do_ad = "true" if config.get("ad_filter", {}).get("enabled", True) else "false",
        ad_expr = _ad_expr(),
        split_bi = "true" if config.get("ad_filter", {}).get("split_multiallelic", True) else "false",
        drop_hets = "true" if config.get("ad_filter", {}).get("drop_hets", True) else "false",
        sample_max_missing = float(config.get("final", {}).get("sample_max_missing", 0.2)),
        site_max_missing = float(config.get("final", {}).get("max_missing", 0.5))
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FINAL} {LOGS}/final

        # help bcftools find plugins in conda envs
        export BCFTOOLS_PLUGINS="$(dirname $(dirname $(which bcftools)))/libexec/bcftools"

        # Step 1: keep PASS and unfiltered (.)
        bcftools view -f .,PASS {input.vcf} -Oz -o {output.passvcf} 2>> {log}
        tabix -p vcf {output.passvcf} 2>> {log}

        # Short-circuit if empty after Step 1
        if [ "$(bcftools view -H {output.passvcf} | head -n1 | wc -l)" -eq 0 ]; then
          echo "[INFO] No SNPs after Step 1; writing empty final VCF." >> {log}
          bcftools view -h {output.passvcf} -Oz -o {output.vcf} 2>> {log}
          tabix -p vcf {output.vcf} 2>> {log}
          cp {output.vcf} {output.advcf}; tabix -p vcf {output.advcf} 2>> {log} || true
          cp {output.vcf} {output.step3vcf}; tabix -p vcf {output.step3vcf} 2>> {log} || true
          : > {output.removes}
          exit 0
        fi

        # Step 2: AD-masking (optionally split multiallelic → biallelic)
        SRC1={output.passvcf}
        if [ "{params.split_bi}" = "true" ]; then
          bcftools norm -m -both "$SRC1" -Oz -o {output.passvcf}.biallelic.gz 2>> {log}
          tabix -p vcf {output.passvcf}.biallelic.gz 2>> {log}
          SRC1={output.passvcf}.biallelic.gz
        fi

        if [ "{params.do_ad}" = "true" ]; then
          if ! bcftools +setGT "$SRC1" -Oz -o {output.advcf} -- -n . -i '{params.ad_expr}' 2>> {log}; then
            echo "[WARN] +setGT (AD mask) failed; skipping AD-clean." >> {log}
            cp "$SRC1" {output.advcf}
          fi
        else
          cp "$SRC1" {output.advcf}
        fi

        # Optionally drop all heterozygous GTs (useful for haploid data)
        if [ "{params.drop_hets}" = "true" ]; then
          if ! bcftools +setGT {output.advcf} -Oz -o {output.advcf}.tmp -- -n . -i 'GT="het"' 2>> {log}; then
            echo "[WARN] +setGT (drop hets) failed; keeping previous genotypes." >> {log}
          else
            mv {output.advcf}.tmp {output.advcf}
          fi
        fi

        tabix -p vcf {output.advcf} 2>> {log}
        SRC="{output.advcf}"

        # Step 3: remove samples with > sample_max_missing missingness
        bcftools query -l "$SRC" > {FINAL}/snps.samples.list 2>> {log}
        ns=$(wc -l < {FINAL}/snps.samples.list)
        vcftools --gzvcf "$SRC" --missing-indv --out {FINAL}/snps 2>> {log}
        LC_NUMERIC=C awk -v thr="{params.sample_max_missing}" 'NR>1 && $5>thr {{print $1}}' \
          {FINAL}/snps.imiss > {output.removes}
        [ -s {output.removes} ] || : > {output.removes}
        nr=$(wc -l < {output.removes})

        if [ "$ns" -eq 0 ] || [ "$nr" -ge "$ns" ]; then
          echo "[WARN] Sample removal would drop all ($nr/$ns). Skipping removal." >> {log}
          cp "$SRC" {output.step3vcf}
          tabix -p vcf {output.step3vcf} 2>> {log} || true
        else
          vcftools --gzvcf "$SRC" \
                   --remove {output.removes} \
                   --recode --recode-INFO-all --stdout 2>> {log} \
            | bgzip -c > {output.step3vcf}
          tabix -p vcf {output.step3vcf} 2>> {log}
        fi

        # Step 4: per-site missingness
        vcftools --gzvcf {output.step3vcf} \
                 --max-missing {params.site_max_missing} \
                 --recode --recode-INFO-all --stdout 2>> {log} \
          | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} 2>> {log}
        """

# ---------------- INDELs ----------------
rule final_filter_indels:
    input:
        vcf = f"{INDELS_DIR}/filtered_indels.merged.vcf.gz",
        tbi = f"{INDELS_DIR}/filtered_indels.merged.vcf.gz.tbi"
    output:
        passvcf = temp(f"{FINAL}/indels.step1_filtered.vcf.gz"),
        passtbi = temp(f"{FINAL}/indels.step1_filtered.vcf.gz.tbi"),
        advcf   = temp(f"{FINAL}/indels.step2_adclean.vcf.gz"),
        adtbi   = temp(f"{FINAL}/indels.step2_adclean.vcf.gz.tbi"),
        step3vcf = temp(f"{FINAL}/indels.step3_filtered.vcf.gz"),
        step3tbi = temp(f"{FINAL}/indels.step3_filtered.vcf.gz.tbi"),
        removes = temp(f"{FINAL}/indels.remove_samples.txt"),
        vcf = f"{FINAL}/filtered_indels.final.vcf.gz",
        tbi = f"{FINAL}/filtered_indels.final.vcf.gz.tbi"
    conda: f"{ENVS}/vcfops.yaml"
    threads: 4
    log: f"{LOGS}/final/indels.log"
    params:
        do_ad = "true" if config.get("ad_filter", {}).get("enabled", True) else "false",
        ad_expr = _ad_expr(),
        split_bi = "true" if config.get("ad_filter", {}).get("split_multiallelic", True) else "false",
        drop_hets = "true" if config.get("ad_filter", {}).get("drop_hets", True) else "false",
        sample_max_missing = float(config.get("final", {}).get("sample_max_missing", 0.2)),
        site_max_missing = float(config.get("final", {}).get("max_missing", 0.5))
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FINAL} {LOGS}/final
        export BCFTOOLS_PLUGINS="$(dirname $(dirname $(which bcftools)))/libexec/bcftools"

        # Step 1: keep PASS and .
        bcftools view -f .,PASS {input.vcf} -Oz -o {output.passvcf} 2>> {log}
        tabix -p vcf {output.passvcf} 2>> {log}

        # Short-circuit if empty after Step 1
        if [ "$(bcftools view -H {output.passvcf} | head -n1 | wc -l)" -eq 0 ]; then
          echo "[INFO] No INDELs after Step 1; writing empty final VCF." >> {log}
          bcftools view -h {output.passvcf} -Oz -o {output.vcf} 2>> {log}
          tabix -p vcf {output.vcf} 2>> {log}
          cp {output.vcf} {output.advcf}; tabix -p vcf {output.advcf} 2>> {log} || true
          cp {output.vcf} {output.step3vcf}; tabix -p vcf {output.step3vcf} 2>> {log} || true
          : > {output.removes}
          exit 0
        fi

        # Step 2: AD-masking (optionally split multiallelic → biallelic)
        SRC1={output.passvcf}
        if [ "{params.split_bi}" = "true" ]; then
          bcftools norm -m -both "$SRC1" -Oz -o {output.passvcf}.biallelic.gz 2>> {log}
          tabix -p vcf {output.passvcf}.biallelic.gz 2>> {log}
          SRC1={output.passvcf}.biallelic.gz
        fi

        if [ "{params.do_ad}" = "true" ]; then
          if ! bcftools +setGT "$SRC1" -Oz -o {output.advcf} -- -n . -i '{params.ad_expr}' 2>> {log}; then
            echo "[WARN] +setGT (AD mask) failed; skipping AD-clean." >> {log}
            cp "$SRC1" {output.advcf}
          fi
        else
          cp "$SRC1" {output.advcf}
        fi

        # Optionally drop all heterozygous GTs
        if [ "{params.drop_hets}" = "true" ]; then
          if ! bcftools +setGT {output.advcf} -Oz -o {output.advcf}.tmp -- -n . -i 'GT="het"' 2>> {log}; then
            echo "[WARN] +setGT (drop hets) failed; keeping previous genotypes." >> {log}
          else
            mv {output.advcf}.tmp {output.advcf}
          fi
        fi

        tabix -p vcf {output.advcf} 2>> {log}
        SRC="{output.advcf}"

        # Step 3: remove samples with > sample_max_missing missingness
        bcftools query -l "$SRC" > {FINAL}/indels.samples.list 2>> {log}
        ns=$(wc -l < {FINAL}/indels.samples.list)
        vcftools --gzvcf "$SRC" --missing-indv --out {FINAL}/indels 2>> {log}
        LC_NUMERIC=C awk -v thr="{params.sample_max_missing}" 'NR>1 && $5>thr {{print $1}}' \
          {FINAL}/indels.imiss > {output.removes}
        [ -s {output.removes} ] || : > {output.removes}
        nr=$(wc -l < {output.removes})

        if [ "$ns" -eq 0 ] || [ "$nr" -ge "$ns" ]; then
          echo "[WARN] Sample removal would drop all ($nr/$ns). Skipping removal." >> {log}
          cp "$SRC" {output.step3vcf}
          tabix -p vcf {output.step3vcf} 2>> {log} || true
        else
          vcftools --gzvcf "$SRC" \
                   --remove {output.removes} \
                   --recode --recode-INFO-all --stdout 2>> {log} \
            | bgzip -c > {output.step3vcf}
          tabix -p vcf {output.step3vcf} 2>> {log}
        fi

        # Step 4: per-site missingness
        vcftools --gzvcf {output.step3vcf} \
                 --max-missing {params.site_max_missing} \
                 --recode --recode-INFO-all --stdout 2>> {log} \
          | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} 2>> {log}
        """
